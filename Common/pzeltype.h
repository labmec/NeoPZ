#ifndef PZELTYPEH
#define PZELTYPEH

// $Id: pzeltype.h,v 1.4 2003-11-07 00:41:50 phil Exp $
/**
   @enum MElementType
 * Defines the element types
 * @param EPoint            element 0D - type point        -  associated index 0
 * @param EOned             element 1D - type oned         -  associated index 1
 * @param ETriangle         element 2D - type triangle     -  associated index 2
 * @param EQuadrilateral    element 2D - type quad         -  associated index 3
 * @param ETetraedro        element 3D - type tetraedro    -  associated index 4
 * @param EPiramide         element 3D - type piramide     -  associated index 5
 * @param EPrisma           element 3D - type prisma       -  associated index 6
 * @param ECube             element 3D - type cube         -  associated index 7
 * @param EPolygonal        element ?? - type ??           -  associated index 8
 * @param EInterface        element nD - type interface    -  associated index 9
 * @param EInterfacePoint   element 1D - type interface    -  associated index 10
 * @param EInterfaceLinear  element 2D - type interface    -  associated index 11
 * @param EInterfaceSurface element 3D - type interface    -  associated index 12
 * @param ESubstructure     element nD - type submesh      -  associated index 13
 * @param EGlobLoc          element nD - type global-local -  associated index 14
 * @param EDiscontinuous    element nD - type discontinuous-  associated index 15
 * @param EInterfaceDisc    element nD - type interface    -  associated index 16
 * @param EAgglomerate      element nD - type agglomerate  -  associated index 17
 * @param ENoType           element 0D - type none         -  associated index 18
 */
//                     0      1        2            3             4            5        
enum MElementType { EPoint, EOned, ETriangle, EQuadrilateral, ETetraedro, EPiramide,
//		       6        7         8           9             10
		    EPrisma,  ECube, EPolygonal, EInterface, EInterfacePoint,
//                        11                 12                13          14       
		    EInterfaceLinear, EInterfaceSurface, ESubstructure, EGlobLoc,
//                       15              16             17          18
		    EDiscontinuous, EAgglomerate, ENoType, EInterfaceDisc = EInterface};
                     
#endif
