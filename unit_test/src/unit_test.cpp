#include "mdimage/include/mdimage.hpp"

int main() {

    std::cout << "hello from mdimage! \n";
    
    auto 
        do_2d     = true,
        do_3d     = true,
        do_4d     = true,
		do_matrix = true;

    float 
        v[2] = { 10, 20 };
    vec2 <float> 
        v2_1 (v);

    if ( do_2d ) {      
        
        std::cout << "2D\nvector v2_1: \t";
        std::cout << v2_1                                           << "\n";
        
        auto v2_2 ( v2_1 );
        v2_2 += 1.;

        std::cout << "vector v2_2: \t"   << v2_2                    << " \n";
        std::cout << "v2_2 + 2 v2_1: \t" << v2_2 + v2_1 * float(2.) << " \n";

    }

    if( do_3d ) {

		vec3 <float>
			v3_1 (v2_1);
        v3_1 [2] = 30;

        std::cout << "3D\nvector v3_1: \t";
        std::cout << v3_1                                           << " \n";

        auto v3_2 ( v3_1 );
        v3_2 += 1.;
        std::cout << "vector v3_2: \t"   << v3_2                    << " \n";
        std::cout << "v3_2 + 2 v3_1: \t" << v3_2 + v3_1 * float(2.) << " \n";

		if ( do_4d ) {

			vec4 <float>
				v4_1 (v3_1);
			v4_1 [3] = 40;

			std::cout << "4D\nvector v4_1: \t";
			std::cout << v4_1                                           << " \n";

			auto v4_2 ( v4_1 );
			v4_2 += 1.;
			std::cout << "vector v4_2: \t"   << v4_2                    << " \n";
			std::cout << "v4_2 + 2 v4_1: \t" << v4_2 + v4_1 * float(2.) << " \n";

		} // if do_4d

    } // if do_3d

    return 0;

}
