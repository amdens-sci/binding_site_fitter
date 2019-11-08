#include <stdio.h>
#include <iostream>
#include <complex.h>
#include "libs\eigen\Eigen\Dense"


extern "C" void find_roots( const double *input_arr, int num_rows, int num_col, double *complex_output_array );

extern "C"
{
  void find_roots( const double *input_array, int num_rows, int num_col,
						double *complex_output_array )
  {
    int num_roots = num_col - 1;
    Eigen::MatrixXd companion_polymat(num_roots, num_roots);
    Eigen::VectorXcd eigenvalues(num_roots);
    Eigen::VectorXd real_component(num_roots);
    Eigen::VectorXd imaginary_component(num_roots);
    for (int input_row=0; input_row < num_rows; ++input_row)
    {
      for (int polymat_row=0; polymat_row < num_roots; ++polymat_row)
      {
        for (int polymat_col=0; polymat_col < num_roots; ++polymat_col)
        {
          if(polymat_row == polymat_col + 1)
            companion_polymat(polymat_row, polymat_col) = 1.0;
          if(polymat_col == num_roots-1)
            companion_polymat(polymat_row, polymat_col) = -input_array[input_row*num_col + polymat_row] / input_array[input_row*num_col + num_roots];
        }
      }
      eigenvalues = companion_polymat.eigenvalues();

      real_component = eigenvalues.real();
      imaginary_component = eigenvalues.imag();
        for(int j=0; j < num_roots; ++j)
        {
          complex_output_array[2*input_row*num_roots + 2*j] = real_component(j);
          complex_output_array[2*input_row*num_roots + 2*j + 1] = imaginary_component(j);
        }
    }
  }


}