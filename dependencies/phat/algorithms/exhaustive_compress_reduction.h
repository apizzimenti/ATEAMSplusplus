/*  Copyright 2013 IST Austria
    Contributed by: Ulrich Bauer, Talha Bin Masood, Barbara Giunti, Guillaume Houry, Michael Kerber, Abhishek Rathod

    This file is part of PHAT.

    PHAT is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PHAT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with PHAT.  If not, see <http://www.gnu.org/licenses/>. */

#pragma once

#include <phat/helpers/misc.h>
#include <phat/boundary_matrix.h>

namespace phat {
  
    class exhaustive_compress_reduction {
	
    public:
	template <typename Representation>
        void operator()(boundary_matrix<Representation> &boundary_matrix) {
	    const index nr_columns = boundary_matrix.get_num_cols();
            std::vector<index> lowest_one_lookup(nr_columns, -1);
	    
	    std::vector<bool> negative(nr_columns,false);
	    
            for (index cur_col = 0; cur_col < nr_columns; cur_col++) {
	        // Compress
	        column faces;
                boundary_matrix.get_col(cur_col, faces);
                column boundary;
                for( index face : faces ) {
                    if(!negative[face]) {
                        boundary.push_back(face);
                    }
                }
                boundary_matrix.set_col(cur_col, boundary);
		
                column final_col_value;
                bool is_pivot = false;
                index lowest_one = boundary_matrix.get_max_index(cur_col);
                while (lowest_one != -1) {
                    if (lowest_one_lookup[lowest_one] != -1) {
                        boundary_matrix.add_to(lowest_one_lookup[lowest_one], cur_col);
                    } else {
                        if (!is_pivot) {
                            lowest_one_lookup[lowest_one] = cur_col;
                            is_pivot = true;
			    negative[cur_col]=true;
                        }
                        final_col_value.push_back(lowest_one);
                        boundary_matrix.remove_max(cur_col);
                    }
                    lowest_one = boundary_matrix.get_max_index(cur_col);
                }
		
                column final_col_vector(final_col_value.rbegin(), final_col_value.rend());
		
                boundary_matrix.set_col(cur_col, final_col_vector);
                boundary_matrix.finalize(cur_col);
            }
        }
    };
}
