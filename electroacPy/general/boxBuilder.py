"""
Functions to automate the creation of simple shoebox-like geometries
"""

import gmsh


class shoebox:
    def __init__(self, Lx, Ly, Lz, position="center", 
                 minSize=0.0057, maxSize=0.057):
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
        self.Lx = Lx
        self.Ly = Ly
        self.Lz = Lz
        self.position = position
        self.minSize = minSize
        self.maxSize = maxSize
        self.membrane = []
        self.name = []


    def addCircularMembrane(self, face, x, y, radius, physical_group, name=None):
        self.membrane.append([face, x, y, radius, physical_group])
        self.name.append(name)
              
    def addRectangularMembrane(self, face, x, y, lx, ly, physical_group, name=None):
        self.membrane.append([face, x, y, lx, ly, physical_group])
        self.name.append(name)
        
    def addPolygon(self, face, X, Y, physical_group, name=None):
        self.membrane.append([face, X, Y, physical_group])
        self.name.append(name)
        
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
        
        ## BUILD PANELS
        for i in range(len(args)):
            piston_tmp = args[i]
            if len(piston_tmp) == 5:
                piston_type = "circle"
            elif len(piston_tmp) == 6:
                piston_type = "rectangle"
            elif len(piston_tmp) == 4:
                piston_type = "polygon"
            
            if piston_type == "circle":
                face = piston_tmp[0]
                x    = piston_tmp[1]
                y    = piston_tmp[2]
                r    = piston_tmp[3]
        
                if face in ["x", "+x", "X", "+X", "-x", "-X"]:
                    if '-' in face:
                        sign = -1
                        Lx_offset = 0
                        Ly_offset = Ly
                    else:
                        sign = 1
                        Lx_offset = Lx
                        Ly_offset = 0
                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(sign * Lx/2, sign*x, y)
                        point_p2 = gmsh.model.geo.addPoint(sign * Lx/2, sign*(x+r), y)
                        point_p3 = gmsh.model.geo.addPoint(sign * Lx/2, sign*x, y+r)
                        point_p4 = gmsh.model.geo.addPoint(sign * Lx/2, sign*(x-r), y)
                        point_p5 = gmsh.model.geo.addPoint(sign * Lx/2, sign*x, y-r)
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lx_offset, Ly_offset + sign*x, y)
                        point_p2 = gmsh.model.geo.addPoint(Lx_offset, Ly_offset + sign*(x+r), y)
                        point_p3 = gmsh.model.geo.addPoint(Lx_offset, Ly_offset + sign*x, y+r)
                        point_p4 = gmsh.model.geo.addPoint(Lx_offset, Ly_offset + sign*(x-r), y)
                        point_p5 = gmsh.model.geo.addPoint(Lx_offset, Ly_offset + sign*x, y-r)
                    circle_1 = gmsh.model.geo.addCircleArc(point_p2, point_p1, point_p3)
                    circle_2 = gmsh.model.geo.addCircleArc(point_p3, point_p1, point_p4)
                    circle_3 = gmsh.model.geo.addCircleArc(point_p4, point_p1, point_p5)
                    circle_4 = gmsh.model.geo.addCircleArc(point_p5, point_p1, point_p2)
                    piston_loop = gmsh.model.geo.addCurveLoop([circle_1, circle_2, 
                                                               circle_3, circle_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                    
                    if "-" in face:
                        panel_xm.append(piston_loop)
                    else:
                        panel_xp.append(piston_loop)
                            
                elif face in ["y", "+y", "Y", "+Y", "-y", "-Y"]:
                    if "-" in face:
                        sign = -1
                        Lx_offset = 0
                        Ly_offset = 0
                    else:
                        sign = 1
                        Lx_offset = Lx
                        Ly_offset = Ly
                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(-sign*x, sign*Ly/2, y)
                        point_p2 = gmsh.model.geo.addPoint(-sign*(x+r), sign*Ly/2, y)
                        point_p3 = gmsh.model.geo.addPoint(-sign*x, sign*Ly/2, y+r)
                        point_p4 = gmsh.model.geo.addPoint(-sign*(x-r), sign*Ly/2, y)
                        point_p5 = gmsh.model.geo.addPoint(-sign*x, sign*Ly/2, y-r)
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lx_offset-sign*x, Ly_offset, y)
                        point_p2 = gmsh.model.geo.addPoint(Lx_offset-sign*(x+r), Ly_offset, y)
                        point_p3 = gmsh.model.geo.addPoint(Lx_offset-sign*x, Ly_offset, y+r)
                        point_p4 = gmsh.model.geo.addPoint(Lx_offset-sign*(x-r), Ly_offset, y)
                        point_p5 = gmsh.model.geo.addPoint(Lx_offset-sign*x, Ly_offset, y-r)
                    circle_1 = gmsh.model.geo.addCircleArc(point_p2, point_p1, point_p3)
                    circle_2 = gmsh.model.geo.addCircleArc(point_p3, point_p1, point_p4)
                    circle_3 = gmsh.model.geo.addCircleArc(point_p4, point_p1, point_p5)
                    circle_4 = gmsh.model.geo.addCircleArc(point_p5, point_p1, point_p2)
                    piston_loop = gmsh.model.geo.addCurveLoop([circle_1, circle_2, 
                                                               circle_3, circle_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                    
                    if "-" in face:
                        panel_ym.append(piston_loop)         
                    else:
                        panel_yp.append(piston_loop)            
        
        
                elif face in ["z", "+z", "Z", "+Z", "-z", "-Z"]:
                    if "-" in face:
                        sign = -1
                        Lz_offset = 0
                        Lx_offset = Lx
                    else:
                        sign = 1
                        Lz_offset = Lz
                        Lx_offset = 0
                                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(sign*x, y, sign*Lz/2)
                        point_p2 = gmsh.model.geo.addPoint(sign*(x+r), y, sign*Lz/2)
                        point_p3 = gmsh.model.geo.addPoint(sign*x, (y+r), sign*Lz/2)
                        point_p4 = gmsh.model.geo.addPoint(sign*(x-r), y, sign*Lz/2)
                        point_p5 = gmsh.model.geo.addPoint(sign*x, (y-r), sign*Lz/2)
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lx_offset+sign*x, y, Lz_offset)
                        point_p2 = gmsh.model.geo.addPoint(Lx_offset+sign*(x+r), y, Lz_offset)
                        point_p3 = gmsh.model.geo.addPoint(Lx_offset+sign*x, (y+r), Lz_offset)
                        point_p4 = gmsh.model.geo.addPoint(Lx_offset+sign*(x-r), y, Lz_offset)
                        point_p5 = gmsh.model.geo.addPoint(Lx_offset+sign*x, y-r, Lz_offset)
                    circle_1 = gmsh.model.geo.addCircleArc(point_p2, point_p1, point_p3)
                    circle_2 = gmsh.model.geo.addCircleArc(point_p3, point_p1, point_p4)
                    circle_3 = gmsh.model.geo.addCircleArc(point_p4, point_p1, point_p5)
                    circle_4 = gmsh.model.geo.addCircleArc(point_p5, point_p1, point_p2)
                    piston_loop = gmsh.model.geo.addCurveLoop([circle_1, circle_2, 
                                                               circle_3, circle_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                    
                    if "-" in face:
                        panel_zm.append(piston_loop)         
                    else:
                        panel_zp.append(piston_loop) 
            
        
            elif piston_type == "rectangle":                
                face = piston_tmp[0]
                x    = piston_tmp[1]
                y    = piston_tmp[2]
                lx   = piston_tmp[3]
                ly   = piston_tmp[4]                
                
                if face in ["x", "+x", "X", "+X", "-x", "-X"]:
                    if '-'  in face:
                        Lxc_offset = -Lx/2
                        Lxr_offset = 0
                        Lyr_offset = Ly
                        sign = -1
                    else:
                        Lxc_offset = Lx/2
                        Lxr_offset = Lx
                        Lyr_offset = 0
                        sign = 1
                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(Lxc_offset, sign*(x-lx/2), y-ly/2)
                        point_p2 = gmsh.model.geo.addPoint(Lxc_offset, sign*(x+lx/2), y-ly/2)
                        point_p3 = gmsh.model.geo.addPoint(Lxc_offset, sign*(x+lx/2), y+ly/2)
                        point_p4 = gmsh.model.geo.addPoint(Lxc_offset, sign*(x-lx/2), y+ly/2)
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lxr_offset, Lyr_offset + sign*x, y)
                        point_p2 = gmsh.model.geo.addPoint(Lxr_offset, Lyr_offset + sign*(x+lx), y)
                        point_p3 = gmsh.model.geo.addPoint(Lxr_offset, Lyr_offset + sign*(x+lx), y+ly)
                        point_p4 = gmsh.model.geo.addPoint(Lxr_offset, Lyr_offset + sign*x, y+ly)
                    
                    rect_1 = gmsh.model.geo.addLine(point_p1, point_p2)
                    rect_2 = gmsh.model.geo.addLine(point_p2, point_p3)
                    rect_3 = gmsh.model.geo.addLine(point_p3, point_p4)
                    rect_4 = gmsh.model.geo.addLine(point_p4, point_p1)
                    piston_loop = gmsh.model.geo.addCurveLoop([rect_1, rect_2, 
                                                               rect_3, rect_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                    
                    if "-" in face:
                        panel_xm.append(piston_loop)         
                    else:
                        panel_xp.append(piston_loop)
        
                
                if face in ["y", "+y", "Y", "+Y", "-y", "-Y"]:
                    if '-' in face:
                        Lyc_offset = -Ly/2
                        Lyr_offset = 0
                        Lxr_offset = 0
                        sign = -1
                    else:
                        Lyc_offset = Ly/2
                        Lyr_offset = Ly
                        Lxr_offset = Lx
                        sign = 1
                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(sign*(x+lx/2), Lyc_offset, y-ly/2)
                        point_p2 = gmsh.model.geo.addPoint(sign*(x-lx/2), Lyc_offset, y-ly/2)
                        point_p3 = gmsh.model.geo.addPoint(sign*(x-lx/2), Lyc_offset, y+ly/2)
                        point_p4 = gmsh.model.geo.addPoint(sign*(x+lx/2), Lyc_offset, y+ly/2)
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lxr_offset - sign*x, Lyr_offset, y)
                        point_p2 = gmsh.model.geo.addPoint(Lxr_offset - sign*(x+lx), Lyr_offset, y)
                        point_p3 = gmsh.model.geo.addPoint(Lxr_offset - sign*(x+lx), Lyr_offset, y+ly)
                        point_p4 = gmsh.model.geo.addPoint(Lxr_offset - sign*x, Lyr_offset, y+ly)
                    
                    rect_1 = gmsh.model.geo.addLine(point_p1, point_p2)
                    rect_2 = gmsh.model.geo.addLine(point_p2, point_p3)
                    rect_3 = gmsh.model.geo.addLine(point_p3, point_p4)
                    rect_4 = gmsh.model.geo.addLine(point_p4, point_p1)
                    piston_loop = gmsh.model.geo.addCurveLoop([rect_1, rect_2, 
                                                               rect_3, rect_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                     
                    if "-" in face:
                        panel_ym.append(piston_loop)         
                    else:
                        panel_yp.append(piston_loop)   
                        
                if face in ["z", "+z", "Z", "+Z", "-z", "-Z"]:                    
                    if '-' in face:
                        Lzc_offset = -Lz/2
                        Lzr_offset = 0
                        Lxr_offset = Lx
                        sign = -1 
                    else:
                        Lzc_offset = Lz/2
                        Lzr_offset = Lz
                        Lxr_offset = 0
                        sign = 1
                        
                    if self.position == "center":
                        point_p1 = gmsh.model.geo.addPoint(sign*(x+lx/2), y-ly/2, Lzc_offset)
                        point_p2 = gmsh.model.geo.addPoint(sign*(x-lx/2), y-ly/2, Lzc_offset)
                        point_p3 = gmsh.model.geo.addPoint(sign*(x-lx/2), y+ly/2, Lzc_offset)
                        point_p4 = gmsh.model.geo.addPoint(sign*(x+lx/2), y+ly/2, Lzc_offset)
                        pos_sign = -1
                    elif self.position == "corner":
                        point_p1 = gmsh.model.geo.addPoint(Lxr_offset + sign*x, y, Lzr_offset)
                        point_p2 = gmsh.model.geo.addPoint(Lxr_offset + sign*(x+lx), y, Lzr_offset)
                        point_p3 = gmsh.model.geo.addPoint(Lxr_offset + sign*(x+lx), y+ly, Lzr_offset)
                        point_p4 = gmsh.model.geo.addPoint(Lxr_offset + sign*x, y+ly, Lzr_offset)
                        pos_sign = 1
                    
                    rect_1 = gmsh.model.geo.addLine(point_p1, point_p2)
                    rect_2 = gmsh.model.geo.addLine(point_p2, point_p3)
                    rect_3 = gmsh.model.geo.addLine(point_p3, point_p4)
                    rect_4 = gmsh.model.geo.addLine(point_p4, point_p1)
                    piston_loop = gmsh.model.geo.addCurveLoop([pos_sign*rect_1, pos_sign*rect_2, 
                                                               pos_sign*rect_3, pos_sign*rect_4])
                    radSurf.append(gmsh.model.geo.addPlaneSurface([piston_loop]))
                     
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
        
        ## SURFACE GROUP
        for i in range(len(args)):
            piston_tmp = args[i]
            if self.name[i] is None:
                name = f"rad_surf_{i+1}"
            else:
                name = self.name[i]
            gmsh.model.addPhysicalGroup(2, [radSurf[i]], tag=piston_tmp[-1], name=name)
        
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
        
        