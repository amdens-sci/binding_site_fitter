#include <stdio.h>
#include <iostream>
#include <complex.h>
#include "../eigen/Eigen/Dense"
#include <cmath>


extern "C" void find_roots( const double *input_arrary, int num_rows, int num_col, double *complex_output_array );
extern "C" void find_single_positive_root( const double *input_arrary, int num_rows, int num_col, double *complex_output_array );

extern "C"
{
  void find_roots( const double *input_array, int num_rows, int num_col,
						double *complex_output_array )
  {
    int num_roots = num_col - 1;
    Eigen::MatrixXd companion_polymat(num_roots, num_roots);
    Eigen::VectorXcd eigenvalues(num_roots);
    for (int input_row=0; input_row < num_rows; ++input_row)
    {
      for (int polymat_row=0; polymat_row < num_roots; ++polymat_row)
      {
        for (int polymat_col=0; polymat_col < num_roots; ++polymat_col)
        {
          if(polymat_row == polymat_col + 1)
            companion_polymat(polymat_row, polymat_col) = 1.0;
          else if(polymat_col == num_roots-1)
            companion_polymat(polymat_row, polymat_col) = -input_array[input_row*num_col + polymat_row] / input_array[input_row*num_col + num_roots];
          else
            companion_polymat(polymat_row, polymat_col) = 0.0;
        }
      }
      eigenvalues = companion_polymat.eigenvalues();

        for(int j=0; j < num_roots; ++j)
        {
          complex_output_array[2*input_row*num_roots + 2*j] = eigenvalues[j].real();
          complex_output_array[2*input_row*num_roots + 2*j + 1] = eigenvalues[j].imag();
        }
    }
  }



//IMPORTANT WARNING: This function should ONLY EVER BE USED if the polynomial you are
//passing can be guaranteed to have only a single positive real root. If the polynomial
//you are passing has multiple positive real roots or no positive real roots, this
//function WILL show highly counterintuitive behavior. I wrote this specifically for a use
//case where the polynomial we are passing can have only a single positive real root and
//the lead coefficient of the polynomial is 1.
//Note that because of issues with floating point arithmetic, we cannot assume that
//something is real only if value.imag() == 0, because there could be cases where 
//value.imag() < 1e-16 and essentially the value IS real but will fail the above test.
//So, that's why we've adopted the procedure shown here.
  void single_pos_special_cubic( const double *input_array, int num_rows, int num_col,
						double *simple_output_array )
  {
    const int num_roots = 3;
    Eigen::Matrix<double, num_roots, num_roots> companion_polymat;
    Eigen::Vector3cd eigenvalues;
    int best_match = 0;
    double smallest_imaginary_value_so_far;
    for (int input_row=0; input_row < num_rows; ++input_row)
    {
      for (int polymat_row=0; polymat_row < num_roots; ++polymat_row)
      {
        for (int polymat_col=0; polymat_col < num_roots; ++polymat_col)
        {
          if(polymat_row == polymat_col + 1)
            companion_polymat(polymat_row, polymat_col) = 1.0;
          else if(polymat_col == num_roots-1)
            companion_polymat(polymat_row, polymat_col) = -input_array[input_row*num_col + polymat_row];
          else
            companion_polymat(polymat_row, polymat_col) = 0.0;
        }
      }
      eigenvalues = companion_polymat.eigenvalues();
      smallest_imaginary_value_so_far = 1;
        for(int j=0; j < num_roots; ++j)
        {
          if( eigenvalues[j].real() > 0 ){
            if( abs(eigenvalues[j].imag()) < smallest_imaginary_value_so_far ){
              smallest_imaginary_value_so_far = abs(eigenvalues[j].imag());
              best_match = j;
            }
          }
        }
      simple_output_array[input_row] = eigenvalues[best_match].real();
    }
  }


//IMPORTANT WARNING: This function should ONLY EVER BE USED if the polynomial you are
//passing can be guaranteed to have only a single positive real root. If the polynomial
//you are passing has multiple positive real roots or no positive real roots, this
//function WILL show highly counterintuitive behavior. I wrote this specifically for a use
//case where the polynomial we are passing can have only a single positive real root and
//the lead coefficient of the polynomial is 1.
//Note that because of issues with floating point arithmetic, we cannot assume that
//something is real only if value.imag() == 0, because there could be cases where 
//value.imag() < 1e-16 and essentially the value IS real but will fail the above test.
//So, that's why we've adopted the procedure shown here.
  void single_pos_special_quartic( const double *input_array, int num_rows, int num_col,
						double *simple_output_array )
  {
    const int num_roots = 4;
    Eigen::Matrix<double, num_roots, num_roots> companion_polymat;
    Eigen::Vector4cd eigenvalues;
    int best_match = 0;
    double smallest_imaginary_value_so_far;
    for (int input_row=0; input_row < num_rows; ++input_row)
    {
      for (int polymat_row=0; polymat_row < num_roots; ++polymat_row)
      {
        for (int polymat_col=0; polymat_col < num_roots; ++polymat_col)
        {
          if(polymat_row == polymat_col + 1)
            companion_polymat(polymat_row, polymat_col) = 1.0;
          else if(polymat_col == num_roots-1)
            companion_polymat(polymat_row, polymat_col) = -input_array[input_row*num_col + polymat_row];
          else
            companion_polymat(polymat_row, polymat_col) = 0.0;
        }
      }
      eigenvalues = companion_polymat.eigenvalues();
      smallest_imaginary_value_so_far = 1;
        for(int j=0; j < num_roots; ++j)
        {
          if( eigenvalues[j].real() > 0 ){
            if( abs(eigenvalues[j].imag()) < smallest_imaginary_value_so_far ){
              smallest_imaginary_value_so_far = abs(eigenvalues[j].imag());
              best_match = j;
            }
          }
        }
      simple_output_array[input_row] = eigenvalues[best_match].real();
    }
  }
}
