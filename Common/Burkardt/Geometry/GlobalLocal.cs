namespace Burkardt.Geometry;

public static class GlobalLocal
{
    public static void glob2loc_3d(double cospitch, double cosroll, double cosyaw,
            double sinpitch, double sinroll, double sinyaw, double[] globas,
            double[] glopts, ref double[] locpts)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GLOB2LOC_3D converts from a global to a local coordinate system in 3D.
        //
        //  Discussion:
        //
        //    A global coordinate system is given.
        //
        //    A local coordinate system has been translated to the point with
        //    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
        //    a roll.
        //
        //    A point has global coordinates GLOPTS, and it is desired to know
        //    the point's local coordinates LOCPTS.
        //
        //    The transformation may be written as
        //
        //      LOC = M_ROLL * M_PITCH * M_YAW * ( GLOB - GLOBAS )
        //
        //    where
        //
        //               (       1            0            0      )
        //    M_ROLL =   (       0        cos(Roll)    sin(Roll)  )
        //               (       0      - sin(Roll)    cos(Roll)  )
        //
        //               (   cos(Pitch)       0      - sin(Pitch) )
        //    M_PITCH =  (       0            1            0      )
        //               (   sin(Pitch)       0        cos(Pitch) )
        //
        //               (   cos(Yaw)     sin(Yaw)         0      )
        //    M_YAW    = ( - sin(Yaw)     cos(Yaw)         0      )
        //               (       0            0            1      )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
        //    roll and yaw angles.
        //
        //    Input, double SINPITCH, SINROLL, SINYAW, the sines of the pitch,
        //    roll and yaw angles.
        //
        //    Input, double GLOBAS[3], the global coordinates of the base vector.
        //
        //    Input, double GLOPTS[3], the global coordinates of the point.
        //
        //    Output, double LOCPTS[3], the local coordinates of the point.
        //
    {
        locpts[0] = cosyaw * cospitch * (glopts[0] - globas[0])
                    + sinyaw * cospitch * (glopts[1] - globas[1])
                    - sinpitch * (glopts[2] - globas[2]);

        locpts[1] = (cosyaw * sinpitch * sinroll - sinyaw * cosroll)
                    * (glopts[0] - globas[0])
                    + (sinyaw * sinpitch * sinroll + cosyaw * cosroll)
                    * (glopts[1] - globas[1])
                    + cospitch * sinroll * (glopts[2] - globas[2]);

        locpts[2] = (cosyaw * sinpitch * cosroll + sinyaw * sinroll)
                    * (glopts[0] - globas[0])
                    + (sinyaw * sinpitch * cosroll - cosyaw * sinroll)
                    * (glopts[1] - globas[1])
                    + cospitch * cosroll * (glopts[2] - globas[2]);

    }

    public static void loc2glob_3d ( double cospitch, double cosroll, double cosyaw,
            double sinpitch, double sinroll, double sinyaw, double[] locpts,
            double[] globas, ref double[] glopts )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOC2GLOB_3D converts from a local to global coordinate system in 3D.
        //
        //  Discussion:
        //
        //    A global coordinate system is given.
        //
        //    A local coordinate system has been translated to the point with
        //    global coordinates GLOBAS, and rotated through a yaw, a pitch, and
        //    a roll.
        //
        //    A point has local coordinates LOCPTS, and it is desired to know
        //    the point's global coordinates GLOPTS.
        //
        //    The transformation may be written as
        //
        //      GLOB = GLOBAS + N_YAW * N_PITCH * N_ROLL * LOC
        //
        //    where
        //
        //               (  cos(Yaw)   -sin(Yaw)        0      )
        //    N_YAW    = (  sin(Yaw)    cos(Yaw)        0      )
        //               (      0           0           1      )
        //
        //               (  cos(Pitch)      0       sin(Pitch) )
        //    N_PITCH =  (      0           1           0      )
        //               ( -sin(Pitch)      0       cos(Pitch) )
        //
        //               (      1           0           0      )
        //    N_ROLL =   (      0       cos(Roll)  -sin(Roll)  )
        //               (      0       sin(Roll)   cos(Roll)  )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double COSPITCH, COSROLL, COSYAW, the cosines of the pitch,
        //    roll and yaw angles.
        //
        //    Input, double SINPITCH, SINROLL, SINYAW, the sines of the pitch,
        //    roll and yaw angles.
        //
        //    Input, double LOCPTS[3], the local coordinates of the point.
        //
        //    Input, double GLOBAS[3], the global coordinates of the base vector.
        //
        //    Output, double GLOPTS[3], the global coordinates of the point.
        //
    {
        glopts[0] = globas[0] + cosyaw * cospitch * locpts[0]
                              + (cosyaw * sinpitch * sinroll - sinyaw * cosroll) * locpts[1]
                              + (cosyaw * sinpitch * cosroll + sinyaw * sinroll) * locpts[2];

        glopts[1] = globas[1] + sinyaw * cospitch * locpts[0]
                              + (sinyaw * sinpitch * sinroll + cosyaw * cosroll) * locpts[1]
                              + (sinyaw * sinpitch * cosroll - cosyaw * sinroll) * locpts[2];

        glopts[2] = globas[2] + -sinpitch * locpts[0]
                              + cospitch * sinroll * locpts[1]
                              + cospitch * cosroll * locpts[2];

    }
}