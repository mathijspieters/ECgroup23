/******************************************************************************
 *  Compilation:  javac Matrix.java
 *  Execution:    java Matrix
 *
 *  A bare-bones collection of static methods for manipulating
 *  matrices.
 *
 ******************************************************************************/
import java.util.Arrays;


public class MatrixFunctions {

    // return n-by-n identity matrix I
    public static double[][] identity(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++){
            a[i][i] = 1;
        }
        return a;
    }

    // return x^T y
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0;
        for (int i = 0; i < x.length; i++){
            sum += x[i] * y[i];
        }
        return sum;
    }

    // return B = A^T
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                b[j][i] = a[i][j];
            }
        }
        return b;
    }

    public static void printMatrix(double[][] a){
        int m = a.length;
        int n = a[0].length;
        System.out.println("--------");
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                System.out.print(a[i][j]);
            }
            System.out.print("\n");
        }
    }

    // return c = a + b
    public static double[][] add(double[][] a, double[][] b) {
        // Adds two matrices
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                c[i][j] = a[i][j] + b[i][j];
            }
        }
        return c;
    }

    public static double[] add(double[] a, double[] b){
        // Adds to vectors
        int m = a.length;
        double[] c = new double[m];
        for (int i = 0; i < m; i++){
            c[i] = a[i] + b[i];
        }
        return c;
    }

    // return c = a - b
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                c[i][j] = a[i][j] - b[i][j];
            }
        }
        return c;
    }

    public static double[] multiply(double[] a, double b){
        // Multiplies vector with constant
        int m = a.length;
        double[] c = new double[m];
        for (int i = 0; i < m; i++){
            c[i] = b * a[i];
        }
        return c;
    }

    public static double[][] multiply(double[] a, double[] b){
        int n = a.length;
        double[][] c = new double[n][n];
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                c[i][j] = a[i]*b[j];
            }
        }
        return c;
    }

    // return c = a * b
    public static double[][] multiply(double[][] a, double b) {
        // Adds two matrices
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                c[i][j] = a[i][j] * b;
            }
        }
        return c;
    }

    // return c = a * b
    public static double[][] multiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++){
            for (int j = 0; j < n2; j++){
                for (int k = 0; k < n1; k++){
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }
        return c;
    }

    // matrix-vector multiplication (y = A * x)
    public static double[] multiply(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++){
            for (int j = 0; j < n; j++){
                y[i] += a[i][j] * x[j];
            }
        }
        return y;
    }


    // vector-matrix multiplication (y = x^T A)
    public static double[] multiply(double[] x, double[][] a) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++){
            for (int i = 0; i < m; i++){
                y[j] += a[i][j] * x[i];
            }
        }
        return y;
    }

    public static boolean isSymmetric(double[][] a) {
        int n = a.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (a[i][j] != a[j][i]){
                    return false;
                }
            }
        }
        return true;
    }

    // is symmetric
    public static boolean isSquare(double[][] a) {
        int n = a.length;
        for (int i = 0; i < n; i++) {
            if (a[i].length != n) return false;
        }
        return true;
    }

    public static double normVec(double[] vector){
        // Returns the norm of a vector

        double norm = 0;
        for (int i = 0; i < vector.length; i++){
            norm += Math.abs(Math.pow(vector[i],2));
        }

        norm = Math.sqrt(norm);
        return norm;
    }


    // return Cholesky factor L of psd matrix A = L L^T
    public static double[][] cholesky(double[][] a) {
        if (!isSquare(a)) {
            throw new RuntimeException("Matrix is not square");
        }
        if (!isSymmetric(a)) {
            throw new RuntimeException("Matrix is not symmetric");
        }

        int n  = a.length;
        double[][] l = new double[n][n];

        for (int i = 0; i < n; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += l[i][k] * l[j][k];
                }
                if (i == j){
                    l[i][i] = Math.sqrt(a[i][i] - sum);
                }
                else{
                     l[i][j] = 1.0 / l[j][j] * (a[i][j] - sum);
                }       
            }
            if (l[i][i] <= 0) {
                throw new RuntimeException("Matrix not positive definite");
            }
        }
        return l;
    }

    public static boolean isPositiveDefinite(double[][] a) {
        if (!isSquare(a)) {
            throw new RuntimeException("Matrix is not square");
        }
        if (!isSymmetric(a)) {
            throw new RuntimeException("Matrix is not symmetric");
        }

        int n  = a.length;
        double[][] l = new double[n][n];

        for (int i = 0; i < n; i++)  {
            for (int j = 0; j <= i; j++) {
                double sum = 0.0;
                for (int k = 0; k < j; k++) {
                    sum += l[i][k] * l[j][k];
                }
                if (i == j){
                    l[i][i] = Math.sqrt(a[i][i] - sum);
                }
                else{
                     l[i][j] = 1.0 / l[j][j] * (a[i][j] - sum);
                }       
            }
            if (l[i][i] <= 0) {
                return false;
            }
        }
        return true;
    }


    public static double[] multivariateGaussian(double[] mean, double[][] covariance, GenUtils genUtil){

        double[][] l = cholesky(covariance);

        //double[][] lTranspose = transpose(l);

        double[] randomGaussian = new double[mean.length];
        for(int i=0; i<mean.length; i++){
            randomGaussian[i] = genUtil.sampleGaussian(0, 1);
        }

        double[] multiVariate = multiply(l,randomGaussian);


        for(int i=0; i<mean.length; i++){
            multiVariate[i] = multiVariate[i] + mean[i];
        }

        return multiVariate;

    }

    public static double[][] matrix_inversion_2(double[][] matrix){
        Matrix m = new Matrix(matrix);
        EigenValueDecomposition e = new EigenValueDecomposition(m);

        double[][] B = e.getV().getArray();
        double[] eigenvalues = e.getRealEigenvalues();

        for(int i=0; i<eigenvalues.length; i++){
            eigenvalues[i] = Math.sqrt(eigenvalues[i]);
            
        }

        double[][] diagonal = new Matrix(eigenvalues).getArray();

        return MatrixFunctions.multiply(B, diagonal);


    }

    public static double[] multivariateGaussianEigenValue(double[] mean, double[][] covariance, GenUtils genUtil){

        double[][] l = matrix_inversion_2(covariance);

        //double[][] lTranspose = transpose(l);

        double[] randomGaussian = new double[mean.length];
        for(int i=0; i<mean.length; i++){
            randomGaussian[i] = genUtil.sampleGaussian(0, 1);
        }

        double[] multiVariate = multiply(l,randomGaussian);


        for(int i=0; i<mean.length; i++){
            multiVariate[i] = multiVariate[i] + mean[i];
        }

        return multiVariate;

    }


}