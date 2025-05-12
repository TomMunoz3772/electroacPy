"""
Functions to automate the creation of simple shoebox-like geometries
"""

import gmsh
import numpy as np


class shoebox:
    """
    Create a parallelepipede.

    Parameters
    ----------
    Lx : float
        Length along x axis.
    Ly : float
        Length along y axis.
    Lz : float
        Length along z axis.
    position : str, optional
        How the box is placed: "center" will put its center at (x, y, z) = 0,
        "corner" will place the lower left corner of the box at (x, y, z) = 0. The default is "center".
    meshSize : float, optional
        Mesh size in metres. The default is 0.057 (343/1000/6).

    Returns
    -------
    None.

    """
    def __init__(self, Lx, Ly, Lz, position="center", 
                 minSize=0.0057, maxSize=0.057):
       
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.position = position
        self.minSize = minSize
        self.maxSize = maxSize
        self.membrane = {}
        # self.name = []
        # self.local = []


    def addCircularBoundary(self, face, x, y, radius, 
                          physical_group, mesh_size=0., name=None, 
                          local="center"):
        """
        Add a circular piston on the given face.

        Parameters
        ----------
        face : str
            Face on which the circular membrane is set. Can be "x", "y", "z", "-x", etc.
        x : float
            Abscissa position of the membrane on the given face.
        y : float
            Ordinate position of the membrane on the given face.
        radius : float
            Radius of the circular membrane.
        physical_group : int
            Physical grouping of the surface. Should be used as reference for the BEM calculations.
        mesh_size : float, optional
            If > 0, size of the membrane's mesh. The default is 0.
        name : str, optional
            Reference given to the physical group. The default is None (will then set 'rad_surf_X').

        Returns
        -------
        None.

        """
        # self.membrane.append([face, x, y, radius, mesh_size, physical_group])
        # self.name.append(name)
        # self.local.append(local)
        if name is None and bool(self.membrane) is False:
            name = "rad_surf_1"
        elif name is None and bool(self.membrane) is True:
            M = len(self.membrane)
            name = f"rad_surf_{M+1}"
        
        self.membrane[name] = circularMembrane(face, x, y, radius,
                                               mesh_size, physical_group, local)

    def addRectangularBoundary(self, face, x, y, lx, ly, 
                             physical_group, mesh_size=0., name=None, 
                             local="center"):
        """
        Add a rectangular membrane on the given face.

        Parameters
        ----------
        face : str
            Face on which the circular membrane is set. Can be "x", "y", "z", "-x", etc.
        x : float
            Abscissa position of the membrane on the given face.
        y : float
            Ordinate position of the membrane on the given face.
        lx : float
            Length along abscissa.
        ly : float
            Length along ordinate.
        physical_group : int
            Physical grouping of the surface. Should be used as reference for the BEM calculations.
        mesh_size : float, optional
            If > 0, size of the membrane's mesh. The default is 0.
        name : str, optional
            Reference given to the physical group. The default is None (will then set 'rad_surf_X').

        Returns
        -------
        None.

        """
        # self.membrane.append([face, x, y, lx, ly, mesh_size, physical_group])
        # self.name.append(name)
        
        if name is None and bool(self.membrane) is False:
            name = "rad_surf_1"
        elif name is None and bool(self.membrane) is True:
            M = len(self.membrane)
            name = f"rad_surf_{M+1}"
        self.membrane[name] = rectangularMembrane(face, x, y, lx, ly, 
                                                  mesh_size, physical_group, local)
        
    def addPolygonBoundary(self, face, X, Y, physical_group, 
                         mesh_size=0., name=None, local="center"):
        """
        Add a polygon as a radiating on the given face.

        Parameters
        ----------
        face : str
            Face on which the circular membrane is set. Can be "x", "y", "z", "-x", etc.
        X : list of float
            Abscissae of all polygon points.
        Y : list of float
            Ordinates of all polygon points. Length should be the same than X.
        physical_group : int
            Physical grouping of the surface. Should be used as reference for the BEM calculations.
        mesh_size : float, optional
            If > 0, size of the membrane's mesh. The default is 0.
        name : str, optional
            Reference given to the physical group. The default is None (will then set 'rad_surf_X').

        Returns
        -------
        None.

        """
        # self.membrane.append([face, X, Y, mesh_size, physical_group])
        # self.name.append(name)
        
        if name is None and bool(self.membrane) is False:
            name = "rad_surf_1"
        elif name is None and bool(self.membrane) is True:
            M = len(self.membrane)
            name = f"rad_surf_{M+1}"
        
        self.membrane[name] = polygonMembrane(face, X, Y, mesh_size, 
                                              physical_group, local)
        
        
    def build(self, path=None):
        """
        Build the geometry.

        Raises
        ------
        ValueError
            Raise error if 'position' is not correctly defined ('center' or 'corner').

        Returns
        -------
        None.

        """
        Lx = self.Lx
        Ly = self.Ly
        Lz = self.Lz
        args = self.membrane
        radSurf = []
        
        gmsh.initialize()
        gmsh.option.setNumber("Mesh.CharacteristicLengthMin", self.minSize)
        gmsh.option.setNumber("Mesh.CharacteristicLengthMax", self.maxSize)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)   

        if self.position == "center":
            point_1 = gmsh.model.geo.addPoint(Lx/2, -Ly/2, -Lz/2)
            point_2 = gmsh.model.geo.addPoint(Lx/2, Ly/2, -Lz/2)
            point_3 = gmsh.model.geo.addPoint(Lx/2, Ly/2, Lz/2)
            point_4 = gmsh.model.geo.addPoint(Lx/2, -Ly/2, Lz/2)
            point_5 = gmsh.model.geo.addPoint(-Lx/2, -Ly/2, -Lz/2)
            point_6 = gmsh.model.geo.addPoint(-Lx/2, Ly/2, -Lz/2)
            point_7 = gmsh.model.geo.addPoint(-Lx/2, Ly/2, Lz/2)
            point_8 = gmsh.model.geo.addPoint(-Lx/2, -Ly/2, Lz/2)
        elif self.position == "corner":
            point_1 = gmsh.model.geo.addPoint(Lx, 0, 0)
            point_2 = gmsh.model.geo.addPoint(Lx, Ly, 0)
            point_3 = gmsh.model.geo.addPoint(Lx, Ly, Lz)
            point_4 = gmsh.model.geo.addPoint(Lx, 0, Lz)
            point_5 = gmsh.model.geo.addPoint(0, 0, 0)
            point_6 = gmsh.model.geo.addPoint(0, Ly, 0)
            point_7 = gmsh.model.geo.addPoint(0, Ly, Lz)
            point_8 = gmsh.model.geo.addPoint(0, 0, Lz)    
        else:
            raise ValueError("'position' should either be 'center' or 'corner'.")
        
        
        ## BUILD LINES
        line_1 = gmsh.model.geo.add_line(point_1, point_2)
        line_2 = gmsh.model.geo.add_line(point_2, point_3)
        line_3 = gmsh.model.geo.add_line(point_3, point_4)
        line_4 = gmsh.model.geo.add_line(point_4, point_1)
        
        line_5 = gmsh.model.geo.add_line(point_2, point_6)
        line_6 = gmsh.model.geo.add_line(point_6, point_7)
        line_7 = gmsh.model.geo.add_line(point_7, point_3)
        line_8 = gmsh.model.geo.add_line(point_7, point_8)

        line_9 = gmsh.model.geo.add_line(point_8, point_4)
        line_10 = gmsh.model.geo.add_line(point_5, point_6)
        line_11 = gmsh.model.geo.add_line(point_5, point_1)
        line_12 = gmsh.model.geo.add_line(point_5, point_8)


        ## MAKE CURVE LOOP FROM LINES
        # curve loop
        loop_1 = gmsh.model.geo.addCurveLoop([line_1, line_2, line_3, line_4])      # +x
        loop_2 = gmsh.model.geo.addCurveLoop([line_5, line_6, line_7, -line_2])     # +y
        loop_3 = gmsh.model.geo.addCurveLoop([-line_3, -line_7, line_8, line_9])      # +z
        loop_4 = gmsh.model.geo.addCurveLoop([-line_11, line_10, -line_5, -line_1])   # -z
        loop_5 = gmsh.model.geo.addCurveLoop([line_11, -line_4, -line_9, -line_12])    # -y
        loop_6 = gmsh.model.geo.addCurveLoop([-line_8, -line_6, -line_10, line_12])   # -x
        
        ## LIST OF PANELS TO BUILD 
        panel_xp = [loop_1]
        panel_xm = [loop_6]
        panel_yp = [loop_2]
        panel_ym = [loop_5]
        panel_zp = [loop_3]
        panel_zm = [loop_4]
        
        
        for m in self.membrane:
            m_tmp = self.membrane[m]
            if m_tmp.mtype == "circular":
                face  = m_tmp.face
                ms    = m_tmp.mesh_size
                xp    = m_tmp.get_placement(Lx, Ly, Lz, self.position)
                
                # place points on face
                point_p1 = gmsh.model.geo.addPoint(xp[0, 0], xp[0, 1], xp[0, 2], ms)
                point_p2 = gmsh.model.geo.addPoint(xp[1, 0], xp[1, 1], xp[1, 2], ms)
                point_p3 = gmsh.model.geo.addPoint(xp[2, 0], xp[2, 1], xp[2, 2], ms)
                point_p4 = gmsh.model.geo.addPoint(xp[3, 0], xp[3, 1], xp[3, 2], ms)
                point_p5 = gmsh.model.geo.addPoint(xp[4, 0], xp[4, 1], xp[4, 2], ms)
                
                # add circle
                circle_1 = gmsh.model.geo.addCircleArc(point_p2, point_p1, point_p3)
                circle_2 = gmsh.model.geo.addCircleArc(point_p3, point_p1, point_p4)
                circle_3 = gmsh.model.geo.addCircleArc(point_p4, point_p1, point_p5)
                circle_4 = gmsh.model.geo.addCircleArc(point_p5, point_p1, point_p2)
                piston_loop = gmsh.model.geo.addCurveLoop([circle_1, circle_2, 
                                                           circle_3, circle_4])
                
                # create new boundary 
                radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                
                if face in ["x", "X", "+x", "+X", "-x", "-X"]:
                    if "-" in face:
                        panel_xm.append(piston_loop)
                    else:
                        panel_xp.append(piston_loop)
                elif face in ["y", "Y", "+y", "+Y", "-y", "-Y"]:
                    if "-" in face:
                        panel_ym.append(piston_loop)
                    else:
                        panel_yp.append(piston_loop)   
                elif face in ["z", "Z", "+z", "+Z", "-z", "-Z"]:
                    if "-" in face:
                        panel_zm.append(piston_loop)
                    else:
                        panel_zp.append(piston_loop)
                        
            if m_tmp.mtype == "rectangular":
                face  = m_tmp.face
                ms    = m_tmp.mesh_size
                xp    = m_tmp.get_placement(Lx, Ly, Lz, self.position)
                                
                # place points on face
                point_p1 = gmsh.model.geo.addPoint(xp[0, 0], xp[0, 1], xp[0, 2], ms)
                point_p2 = gmsh.model.geo.addPoint(xp[1, 0], xp[1, 1], xp[1, 2], ms)
                point_p3 = gmsh.model.geo.addPoint(xp[2, 0], xp[2, 1], xp[2, 2], ms)
                point_p4 = gmsh.model.geo.addPoint(xp[3, 0], xp[3, 1], xp[3, 2], ms)
                
                rect_1 = gmsh.model.geo.addLine(point_p1, point_p2)
                rect_2 = gmsh.model.geo.addLine(point_p2, point_p3)
                rect_3 = gmsh.model.geo.addLine(point_p3, point_p4)
                rect_4 = gmsh.model.geo.addLine(point_p4, point_p1)
                
                piston_loop = gmsh.model.geo.addCurveLoop([rect_1, rect_2, 
                                                           rect_3, rect_4])
                
                # create new boundary 
                radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                
                if face in ["x", "X", "+x", "+X", "-x", "-X"]:
                    if "-" in face:
                        panel_xm.append(piston_loop)
                    else:
                        panel_xp.append(piston_loop)
                elif face in ["y", "Y", "+y", "+Y", "-y", "-Y"]:
                    if "-" in face:
                        panel_ym.append(piston_loop)
                    else:
                        panel_yp.append(piston_loop)   
                elif face in ["z", "Z", "+z", "+Z", "-z", "-Z"]:
                    if "-" in face:
                        panel_zm.append(piston_loop)
                    else:
                        panel_zp.append(piston_loop)
                        
                        
            if m_tmp.mtype == "polygon":
                face  = m_tmp.face
                ms    = m_tmp.mesh_size
                xp = m_tmp.get_placement(Lx, Ly, Lz, self.position)
                point_px = []
                poly = []
                
                for i in range(len(xp)):
                    point_px.append(gmsh.model.geo.addPoint(xp[i, 0],
                                                            xp[i, 1], xp[i, 2]))    
                    
                # link points together (lines)
                for j in range(len(point_px)-1):
                    poly.append(gmsh.model.geo.addLine(point_px[j], point_px[j+1]))
                poly.append(gmsh.model.geo.addLine(point_px[-1], point_px[0]))
                piston_loop = gmsh.model.geo.addCurveLoop(poly)
                radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
            
                if face in ["x", "X", "+x", "+X", "-x", "-X"]:
                    if "-" in face:
                        panel_xm.append(piston_loop)
                    else:
                        panel_xp.append(piston_loop)
                elif face in ["y", "Y", "+y", "+Y", "-y", "-Y"]:
                    if "-" in face:
                        panel_ym.append(piston_loop)
                    else:
                        panel_yp.append(piston_loop)   
                elif face in ["z", "Z", "+z", "+Z", "-z", "-Z"]:
                    if "-" in face:
                        panel_zm.append(piston_loop)
                    else:
                        panel_zp.append(piston_loop)
                        
            
        nonRadSurf = []
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_xp))
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_xm))
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_yp))
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_ym))
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_zp))
        nonRadSurf.append(gmsh.model.geo.addPlaneSurface(panel_zm))
        gmsh.model.geo.synchronize()
        
        # create rad/impedance groups
        groups = {}
        gname  = {}
        for i, m in enumerate(self.membrane):
            m_tmp = self.membrane[m]
            group_tmp = str(m_tmp.physical_group)
            if group_tmp not in groups:
                groups[group_tmp] = [radSurf[i]]
                gname[group_tmp] = m
            else:
                groups[group_tmp].append(radSurf[i])
        
        for g in groups:
            surf_tmp = groups[g]
            gmsh.model.addPhysicalGroup(2, surf_tmp, tag=int(g), name=gname[g])
        
        # remaining group (enclosure / room)
        gmsh.model.addPhysicalGroup(2, nonRadSurf, name="rigid_boundary")
        
        ## BUILD GEOMETRY        
        gmsh.model.mesh.generate(dim=2)
        
        if path is None:
            gmsh.write("geo_mesh.msh")
        elif type(path) == str:
            if path[-4::] in [".med", ".msh"]:
                gmsh.write(path)
            else:
                gmsh.finalize()
                raise Exception("Mesh extension not supported. Try '.med' or '.msh'.")
            
        gmsh.finalize()
        return None
        
        
        
#%% SHAPE INFO HOLDER
class circularMembrane:
    def __init__(self, face, x, y, radius, mesh_size, physical_group, local):
        self.face           = face
        self.x              = x
        self.y              = y
        self.radius         = radius
        self.mesh_size      = mesh_size
        self.physical_group = physical_group
        self.local          = local
        self.mtype          = "circular"
        
    def get_placement(self, Lx, Ly, Lz, position):
        # for more lisibility
        x = self.x
        y = self.y
        r = self.radius
        
        
        if self.face in ["x", "+x", "X", "+X", "-x", "-X"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = -Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = Ly
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = 0
                    Lz_offset = 0
            
            # point position
            xp = np.array([[0, x, y], 
                           [0, x+r, y], 
                           [0, x, y+r], 
                           [0, x-r, y],
                           [0, x, y-r]])
            xp[:, 1] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset
        
        if self.face in ["y", "+y", "Y", "+Y", "-y", "-Y"]:
            # check which side of the box
            if "-" in self.face:
                sign = 1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = -Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = -1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = Ly
                    Lz_offset = 0
            # point position
            xp = np.array([[x, 0, y], 
                           [x+r, 0, y], 
                           [x, 0, y+r], 
                           [x-r, 0, y],
                           [x, 0, y-r]])
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset    
                         
            
        if self.face in ["z", "+z", "Z", "+Z", "-z", "-Z"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = -Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = +Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = Lz
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = Lz
            
            # point position
            xp = np.array([[x, y, 0], 
                           [x+r, y, 0], 
                           [x, y+r, 0], 
                           [x-r, y, 0],
                           [x, y-r, 0]])
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset        
        return xp
        
        
class rectangularMembrane:
    def __init__(self, face, x, y, lx, ly, mesh_size, physical_group, local):
        self.face           = face
        self.x              = x
        self.y              = y
        self.lx             = lx
        self.ly             = ly
        self.mesh_size      = mesh_size
        self.physical_group = physical_group
        self.local          = local
        self.mtype          = "rectangular"
        
    def get_placement(self, Lx, Ly, Lz, position):
        # for more lisibility
        x = self.x
        y = self.y
        lx = self.lx
        ly = self.ly
        
        
        if self.face in ["x", "+x", "X", "+X", "-x", "-X"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = -Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = Ly
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = 0
                    Lz_offset = 0
            
            # point position
            xp = np.array([[0, x-lx/2, y-ly/2], 
                           [0, x+lx/2, y-ly/2], 
                           [0, x+lx/2, y+ly/2], 
                           [0, x-lx/2, y+ly/2]])
            xp[:, 1] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset
        
        if self.face in ["y", "+y", "Y", "+Y", "-y", "-Y"]:
            # check which side of the box
            if "-" in self.face:
                sign = 1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = -Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = -1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = Ly
                    Lz_offset = 0
            # point position
            xp = np.array([[x-lx/2, 0, y-ly/2], 
                           [x+lx/2, 0, y-ly/2], 
                           [x+lx/2, 0, y+ly/2], 
                           [x-lx/2, 0, y+ly/2]])
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset    
                         
            
        if self.face in ["z", "+z", "Z", "+Z", "-z", "-Z"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = -Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = Lz
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = Lz
            
            # point position
            xp = np.array([[x-lx/2, y-ly/2, 0], 
                           [x+lx/2, y-ly/2, 0], 
                           [x+lx/2, y+ly/2, 0], 
                           [x-lx/2, y+ly/2, 0]])
            
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset        
        return xp        
    
        
class polygonMembrane:
    def __init__(self, face, X, Y, mesh_size, physical_group, local):
        self.face           = face
        self.X              = X
        self.Y              = Y
        self.mesh_size      = mesh_size
        self.physical_group = physical_group
        self.local          = local
        self.mtype          = "polygon"
        

    def get_placement(self, Lx, Ly, Lz, position):
        # for more lisibility
        X = self.X
        Y = self.Y
        
        
        if self.face in ["x", "+x", "X", "+X", "-x", "-X"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = -Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = Ly
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = 0
                    Lz_offset = 0
            
            # point position
            xp = np.zeros([len(X), 3])
            
            for i in range(len(X)):
                xp[i, 1] = X[i]
                xp[i, 2] = Y[i]
            
            xp[:, 1] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset
        
        if self.face in ["y", "+y", "Y", "+Y", "-y", "-Y"]:
            # check which side of the box
            if "-" in self.face:
                sign = 1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = -Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = -1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "center" and self.local == "corner":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "corner":
                    Lx_offset = Lx
                    Ly_offset = Ly
                    Lz_offset = 0

            # point position
            xp = np.zeros([len(X), 3])
            
            for i in range(len(X)):
                xp[i, 0] = X[i]
                xp[i, 2] = Y[i]
            
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset
            
            
        if self.face in ["z", "+z", "Z", "+Z", "-z", "-Z"]:
            # check which side of the box
            if "-" in self.face:
                sign = -1            
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = -Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = -Ly/2
                    Lz_offset = -Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = 0
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = 0 
            else:
                sign = 1
                # get local coordinates (Lx_offset, Ly_offset, Lz_offset)                    
                if position == "center" and self.local == "center":
                    Lx_offset = 0
                    Ly_offset = 0
                    Lz_offset = Lz/2
                elif position == "center" and self.local == "corner":
                    Lx_offset = -Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = Lz/2
                elif position == "corner" and self.local == "center":
                    Lx_offset = Lx/2
                    Ly_offset = Ly/2
                    Lz_offset = Lz
                elif position == "corner" and self.local == "corner":
                    Lx_offset = 0
                    Ly_offset = Ly
                    Lz_offset = Lz
            
            # point position
            xp = np.zeros([len(X), 3])
            
            for i in range(len(X)):
                xp[i, 0] = X[i]
                xp[i, 1] = Y[i]
            
            xp[:, 0] *= sign
            xp[:, 0] += Lx_offset
            xp[:, 1] += Ly_offset
            xp[:, 2] += Lz_offset  
        return xp        


