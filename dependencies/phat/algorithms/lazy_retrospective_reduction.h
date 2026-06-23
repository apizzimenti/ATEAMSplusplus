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
#include <stack>

namespace phat {
    class lazy_retrospective_reduction {
    public:
        template< typename Representation >
        void operator () ( boundary_matrix< Representation >& boundary_matrix ) {
	    const index nr_columns = boundary_matrix.get_num_cols();
            std::vector< index > positivePair( nr_columns, -1 );
            std::vector< index > negativePair( nr_columns, -1 );
	    
            for( index j = 0; j < nr_columns; j++ ) {
                if (boundary_matrix.is_empty(j)) {
                    continue;
                }
		
                // Compress
                column faces;
                boundary_matrix.get_col(j, faces);
                column boundary;
                for( index face : faces ) {
                    if(negativePair[face] == -1) {
                        boundary.push_back(face);
                    }
                }
                boundary_matrix.set_col(j, boundary);
		
                std::stack <index> s;
                s.push(j); 
                while (!s.empty()) { 
                    index l = s.top();
                    // column faces;
                    boundary_matrix.get_col(l, faces);
                    bool foundPositivePaired = false;
                    for( index face : faces ) {
                        index pair = positivePair[face];
                        if( pair != -1 && pair != l){
                            s.push(pair);
                            foundPositivePaired = true;
                            break;
                        }
                    }
                    if (!foundPositivePaired) {
                        index t = s.top();
                        s.pop();
                        if (!s.empty()){
                            index l = s.top();
                            boundary_matrix.add_to(t, l);
                        }
                    }
                }
                if (!boundary_matrix.is_empty(j)) {
                    index pair = boundary_matrix.get_max_index(j);
                    positivePair[pair] = j;
                    negativePair[j] = pair;
                }
            }
        }
    };
}
