# GitHub Wiki Deployment Guide

This guide explains how to populate the GitHub project wiki with the comprehensive API documentation created in the `Wiki/` directory.

## Understanding GitHub Wikis

GitHub wikis are **separate Git repositories** from your main project repository. The wiki content in the `Wiki/` directory of this repository will **NOT automatically populate** the GitHub wiki when merged. Manual deployment is required.

## Option 1: Manual Wiki Page Creation (Recommended for First-Time Setup)

### Step 1: Access the GitHub Wiki
1. Go to your repository on GitHub: `https://github.com/philstopford/DesignLibs_GPL`
2. Click the "Wiki" tab at the top
3. If no wiki exists, click "Create the first page"

### Step 2: Create the Home Page
1. Click "Create new page" or edit the existing Home page
2. Copy the entire contents of `Wiki/Home.md` into the page editor
3. Set the page title as "Home" 
4. Click "Save Page"

### Step 3: Create Individual Library Pages
For each `.md` file in the `Wiki/` directory:

1. **geoLib-API.md** → Create new page titled "geoLib API"
2. **geoCore-API.md** → Create new page titled "geoCore API"
3. **geoWrangler-API.md** → Create new page titled "geoWrangler API"
4. **shapeEngine-API.md** → Create new page titled "shapeEngine API"
5. **clipper-API.md** → Create new page titled "clipper API"
6. **utility-API.md** → Create new page titled "utility API"
7. **MersenneTwister-API.md** → Create new page titled "MersenneTwister API"
8. **Noise-API.md** → Create new page titled "Noise API"
9. **Eto.VeldridSurface-API.md** → Create new page titled "Eto.VeldridSurface API"

For each page:
- Click "New Page" in the wiki
- Copy the entire contents from the corresponding `.md` file
- Use the page title format shown above
- Click "Save Page"

## Option 2: Command-Line Deployment (Advanced)

### Prerequisites
- Git command line access
- Write permissions to the repository

### Step 1: Clone the Wiki Repository
```bash
# Clone the separate wiki repository
git clone https://github.com/philstopford/DesignLibs_GPL.wiki.git

# Navigate to the wiki directory
cd DesignLibs_GPL.wiki
```

### Step 2: Copy Documentation Files
```bash
# Copy all markdown files from the main repository's Wiki directory
cp /path/to/DesignLibs_GPL/Wiki/*.md .

# Rename files to match GitHub wiki conventions (spaces become dashes)
mv "geoLib-API.md" "geoLib-API.md"
mv "geoCore-API.md" "geoCore-API.md" 
mv "geoWrangler-API.md" "geoWrangler-API.md"
mv "shapeEngine-API.md" "shapeEngine-API.md"
mv "clipper-API.md" "clipper-API.md"
mv "utility-API.md" "utility-API.md"
mv "MersenneTwister-API.md" "MersenneTwister-API.md"
mv "Noise-API.md" "Noise-API.md"
mv "Eto.VeldridSurface-API.md" "Eto-VeldridSurface-API.md"
```

### Step 3: Commit and Push to Wiki
```bash
# Add all new documentation
git add *.md

# Commit the documentation
git commit -m "Add comprehensive API documentation for all DesignLibs_GPL libraries"

# Push to the wiki repository
git push origin master
```

## Option 3: Automated Deployment Script

I can create a deployment script to automate this process:

```bash
#!/bin/bash
# deploy-wiki.sh

REPO_URL="https://github.com/philstopford/DesignLibs_GPL.wiki.git"
WIKI_SOURCE_DIR="Wiki"
TEMP_WIKI_DIR="/tmp/wiki-deployment"

echo "Deploying DesignLibs_GPL wiki documentation..."

# Clone wiki repository
git clone "$REPO_URL" "$TEMP_WIKI_DIR"
cd "$TEMP_WIKI_DIR"

# Copy documentation files
cp "../$WIKI_SOURCE_DIR"/*.md .

# Add and commit changes
git add *.md
git commit -m "Update API documentation from repository"
git push origin master

echo "Wiki deployment complete!"
echo "Visit https://github.com/philstopford/DesignLibs_GPL/wiki to view"

# Cleanup
cd ..
rm -rf "$TEMP_WIKI_DIR"
```

## Important Notes

### File Naming Conventions
- GitHub wiki pages use the filename as the page title
- Hyphens in filenames become spaces in page titles
- The `Home.md` file becomes the wiki home page

### Link Updates Required
After deploying, you may need to update internal links in the documentation:
- Repository links like `geoLib-API.md` should become `geoLib-API` (remove .md extension)
- Cross-references between pages will work automatically

### Content Verification
After deployment:
1. Visit `https://github.com/philstopford/DesignLibs_GPL/wiki`
2. Verify all pages are created correctly
3. Check that links between pages work properly
4. Confirm formatting appears correctly

## Documentation Statistics

The created documentation includes:
- **10 comprehensive API documentation files** (5,529 total lines)
- **Complete coverage** of all major DesignLibs_GPL libraries
- **Detailed API references** with examples and usage patterns
- **Cross-references** showing library interactions
- **Ready-to-use** content formatted for GitHub wiki

## Maintenance

To update the wiki in the future:
1. Make changes to files in the repository's `Wiki/` directory
2. Re-run the deployment process using your chosen method
3. The wiki will be updated with the latest documentation

This documentation system provides professional-quality API reference material for the entire DesignLibs_GPL ecosystem.