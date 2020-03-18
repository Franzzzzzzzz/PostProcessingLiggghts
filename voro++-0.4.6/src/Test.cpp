// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
using namespace voro;

// Set up constants for the container geometry
const double x_min=0,x_max=1;
const double y_min=0,y_max=1;
const double z_min=-0.5,z_max=0.5;

// Set up the number of blocks that the container is divided
// into.
const int n_x=3,n_y=3,n_z=3;

int main() {

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block. Import
	// the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	//container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
	//		false,false,false,8);
	//con.import("pack_six_cube");
	//con.draw_cells_gnuplot("pack_six_cube.gnu");
	//con.draw_cells_pov("pack_six_cube_v.pov");
	//con.draw_particles_pov("pack_six_cube_p.pov");

	// Create a container for polydisperse particles using the same
	// geometry as above. Import the polydisperse test packing and
	// output the Voronoi radical tessellation in gnuplot and POV-Ray
	// formats.
	container_poly cont(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);
    
    
    cont.put(1,0.001,0.001,0,0.1);
    cont.put(2,0.5  ,0.001,0,0.1);
    cont.put(3,0.001,0.5  ,0,0.1);
    cont.put(4,0.5  ,0.5  ,0,0.1);
    cont.put(5,0.25 ,0.25 ,0,0.25);
    //cont.put(5,0.5 ,0   ,0,0.01);
    
    voronoicell a ;
    c_loop_all b(cont); 
    cont.compute_all_cells(); 
    b.start() ; 
    cont.draw_cells_gnuplot("Test.gpl");
    
    int p ; double u, v, w ; 
        cont.find_voronoi_cell(0.6,0.5,0.5,u,v,w,p) ;
    do
    {
        double x, y, z ; 
        cont.compute_cell(a,b) ;
        a.centroid(x,y,z) ; 
        printf("%g %g %g | ", b.x(), b.y(), b.z()) ; 
        printf("%g %g %g %g %d %g %g %g %d\n", a.volume(), x, y, z, a.p, u, v, w, p) ; 
    } 
    while (b.inc()) ;
    
    /*
	conp.import("pack_six_cube_poly");
	conp.draw_cells_gnuplot("pack_six_cube_poly.gnu");
	conp.draw_cells_pov("pack_six_cube_poly_v.pov");
	conp.draw_particles_pov("pack_six_cube_poly_p.pov");*/
}
