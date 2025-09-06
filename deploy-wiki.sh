#!/bin/bash
# deploy-wiki.sh - Automated GitHub Wiki Deployment Script
# 
# This script deploys the API documentation from the Wiki/ directory 
# to the GitHub wiki for DesignLibs_GPL
#
# Usage: ./deploy-wiki.sh

set -e  # Exit on any error

REPO_NAME="DesignLibs_GPL"
REPO_OWNER="philstopford"
WIKI_URL="https://github.com/${REPO_OWNER}/${REPO_NAME}.wiki.git"
WIKI_SOURCE_DIR="Wiki"
TEMP_WIKI_DIR="/tmp/designlibs-wiki-deployment"

echo "🚀 DesignLibs_GPL Wiki Deployment Script"
echo "========================================="

# Check if Wiki source directory exists
if [ ! -d "$WIKI_SOURCE_DIR" ]; then
    echo "❌ Error: Wiki source directory '$WIKI_SOURCE_DIR' not found!"
    echo "   Please run this script from the repository root directory."
    exit 1
fi

# Count documentation files
WIKI_FILES=$(find "$WIKI_SOURCE_DIR" -name "*.md" | wc -l)
echo "📚 Found $WIKI_FILES documentation files to deploy"

# Clean up any existing temporary directory
if [ -d "$TEMP_WIKI_DIR" ]; then
    echo "🧹 Cleaning up previous deployment directory..."
    rm -rf "$TEMP_WIKI_DIR"
fi

echo "📥 Cloning wiki repository..."
if git clone "$WIKI_URL" "$TEMP_WIKI_DIR" 2>/dev/null; then
    echo "✅ Wiki repository cloned successfully"
else
    echo "❌ Error: Could not clone wiki repository"
    echo "   Make sure you have write access to the repository and the wiki is enabled"
    echo "   You may need to create the first wiki page manually on GitHub first"
    exit 1
fi

# Navigate to wiki directory
cd "$TEMP_WIKI_DIR"

echo "📋 Copying documentation files..."

# Copy all markdown files from the source directory
for file in "../$WIKI_SOURCE_DIR"/*.md; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        cp "$file" "$filename"
        echo "   ✓ Copied $filename"
    fi
done

echo "🔍 Checking for changes..."

# Check if there are any changes to commit
if git diff --quiet && git diff --cached --quiet; then
    echo "ℹ️  No changes detected - wiki is already up to date"
    cd ..
    rm -rf "$TEMP_WIKI_DIR"
    exit 0
fi

# Show what files will be added/modified
echo "📝 Changes to be committed:"
git status --porcelain

# Add all changes
echo "➕ Adding documentation files..."
git add *.md

# Commit with timestamp
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
COMMIT_MSG="Update DesignLibs_GPL API documentation - $TIMESTAMP

- Updated comprehensive API documentation for all major libraries
- Includes geoLib, geoCore, geoWrangler, shapeEngine, clipper, utility, MersenneTwister, Noise, and Eto.VeldridSurface
- Total documentation: $WIKI_FILES files with complete API coverage
- Generated from repository Wiki/ directory"

echo "💾 Committing changes..."
git commit -m "$COMMIT_MSG"

echo "🚀 Pushing to GitHub wiki..."
if git push origin master; then
    echo "✅ Wiki deployment successful!"
    echo ""
    echo "🌐 Your wiki is now available at:"
    echo "   https://github.com/${REPO_OWNER}/${REPO_NAME}/wiki"
    echo ""
    echo "📊 Deployment Summary:"
    echo "   - $WIKI_FILES documentation files deployed"
    echo "   - Complete API coverage for DesignLibs_GPL libraries"
    echo "   - Professional-quality reference documentation"
else
    echo "❌ Error: Could not push changes to wiki repository"
    echo "   Check your permissions and try again"
    cd ..
    rm -rf "$TEMP_WIKI_DIR"
    exit 1
fi

# Clean up
echo "🧹 Cleaning up temporary files..."
cd ..
rm -rf "$TEMP_WIKI_DIR"

echo "🎉 Wiki deployment completed successfully!"
echo "   Visit https://github.com/${REPO_OWNER}/${REPO_NAME}/wiki to view your documentation"