
#include <iostream>
#include <limits.h>
#include <stdio.h>

#define EIGEN2_SUPPORT
#define EIGEN_NO_EIGEN2_DEPRECATED_WARNING
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <eigen3/Eigen/LeastSquares>


using namespace Eigen;
using namespace std;


// project vect onto plane which has e1 and e2 as orthonormal basis
// return result in proj 

static Vector3d project_to_plane(Vector3d &e1, Vector3d &e2, Vector3d &vect)
{
   int row, col;
   Vector3d tmp;

   tmp = ((vect.dot(e1) / e1.dot(e1)) * e1) + (((vect.dot(e2) / e2.dot(e2))) * e2);

//   cout << tmp.rows() << "   " << tmp.cols() << endl;

//   cout << "tiny is: "  << std::numeric_limits<double>::epsilon() << endl;

   for (row = 0 ; row < tmp.rows(); row++)
      for (col = 0 ; col < tmp.cols(); col++)
         if (tmp(row,col) != 0 && std::abs(tmp(row,col)) <= std::numeric_limits<double>::epsilon())
         {
            cout << endl << "change " << tmp(row,col) <<  " tiny to zero " << col << "  " << row << endl << endl;
            cout << "before: " << endl <<  tmp << endl; 
            tmp(row,col) = (double) 0.0;
            cout << "after: " << endl <<  tmp << endl << endl; 
         }

   return tmp;
}

int main(int argc, char** argv)
{
   Vector3d A;
   Vector3d B;
   Vector3d C;
   Vector3d CA;
   Vector3d CB;
   Vector3d norm_cb;
   Vector3d perp_to_cb;
   Vector3d norm_perp_to_cb;
   Vector3d tmp_vect;
   Vector3d proj;
   double  dotter, rhs, lhs;
   double  len_cb;
   double ca_proj_on_cb;

//   cout.precision(8);

   MatrixXd m(3,3);

   m << 1, 2, 3, 2, 4, 6, -1, 0, 12;

   cout << m << endl << endl;

   m.col(0) << 1, 2, 3;
   m.col(1) << 2, 4, 6;
   m.col(2) << -1, 0, 12;
   cout << m << endl << endl;

   A << 4, 7, 0;
   B << 5, 3, 0 ;
   C << 2, 1, 0 ;
   cout << "A:" <<endl << A << endl << endl;
   cout << "B:" <<endl << B << endl << endl;
   cout << "C:" <<endl << C << endl << endl;

      // position vectors
   CB = B - C;
   CA = A - C;

   cout << "CB: " << endl << CB << endl << endl;
   cout << "CA: " << endl << CA << endl << endl;

    // colinear?
   lhs = (CA+CB).norm();
   rhs = CA.norm() + CB.norm();
   cout << "Test for colinear" << endl << "AC:     " << lhs << endl << "AB->AC: " << rhs << endl << endl;
   if (std::abs(lhs - rhs) < std::numeric_limits<double>::epsilon())
   {
      cout << "Points colinear, no plane possible." << endl;
      exit(1);
   }
   else
      cout << "OK!" << endl << endl;


   dotter = CA.dot(CB);

   cout << "dot: " << dotter << endl;
   Vector3d cb_norm = CB;
   Vector3d ca_norm = CA;
   cb_norm.normalize();
   ca_norm.normalize();

   cout << "norm dot: " << ca_norm.dot(cb_norm) << endl << endl;

   len_cb = CB.norm();
   norm_cb = CB / len_cb;

   cout << "CB len  " << len_cb << endl;
   cout << "norm cb" << endl <<  norm_cb << endl << endl;

   ca_proj_on_cb = CA.dot(norm_cb); 
   cout << "ca component proj on cb " << ca_proj_on_cb << endl;

     // find orthonormal basis for plane
   perp_to_cb = CA - (CA.dot(CB) / CB.dot(CB)) * CB;
   cout << "perp to cb" << endl << perp_to_cb << endl << endl;

   norm_perp_to_cb = perp_to_cb / perp_to_cb.norm();
   cout << "norm perp to cb" << endl << norm_perp_to_cb << endl << endl;

   cout << "Check for orth basis:  " << norm_cb.dot(norm_perp_to_cb) << endl << endl;
   if (norm_cb.dot(norm_perp_to_cb) <=  std::numeric_limits<double>::epsilon())
      cout << "Orthoganal" << endl << endl;
   else
      cout << "NOT Orthoganal" << endl << endl;

    // project some vects.  Head will be point projected into plane
   tmp_vect << 5, 6, 7;
   cout << "temp 1: " << endl << tmp_vect << endl;
   proj = project_to_plane(norm_cb, norm_perp_to_cb, tmp_vect);
   cout << "proj 1: " << endl << proj << endl;

   tmp_vect << -4, 0, 8;
   cout << "temp 2: " << endl << tmp_vect << endl;
   proj = project_to_plane(norm_cb, norm_perp_to_cb, tmp_vect);
   cout << "proj 2: " << endl << proj << endl;

   tmp_vect << -11, -34, -97;
   cout << "temp 3: " << endl << tmp_vect << endl;
   proj = project_to_plane(norm_cb, norm_perp_to_cb, tmp_vect);
   cout << "temp : " << endl << proj << endl;


   // block ops

   cout << "head:  " << A.head(1) << endl << endl;



   // solve system

   Matrix3d F;
   Vector3d b;
//   F << 1,2,3,4,5,6,7,8,10;
   F.col(0) = A;
   F.col(1) = B;
   F.col(2) = C;
   b << 0,0,0;

// 4a1 + 7a2 + 0a3 = 0
// 5b1 + 3b2 + 0b3 = 0
// 2a3 + 1y + 0z = 0

   cout << "F:\n" << F << endl;
   cout << "b:\n" << b << endl;

   Vector3d F_res = F.colPivHouseholderQr().solve(b);


   cout << " res is:\n" << F_res << endl << endl;

   F_res = F.fullPivHouseholderQr().solve(b);
   cout << " res #2 is:\n" << F_res << endl << endl;
   double solve;

   double d = F_res(0) * A(0) + F_res(1) * A(1) + F_res(2) * A(2);
   cout << "d: " << d << endl << endl;

   solve = F_res(0) * B(0) + F_res(1) * B(1) + F_res(2) * B(2);
   cout << "d: " << solve << endl << endl;
   return 0;
}


