package beast.base.evolution.tree.coalescent;

import beast.base.core.Description;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.function.Cos;
import org.apache.commons.math3.analysis.function.Sin;

@Description("Creates a truncated Fourier series A0 + sum(An*cos(n*pi*t)+Bn*sin(n*pi*t))")
public class TruncatedFourierSeries implements UnivariateFunction {

    private double[] coefficients;

    private final double period;

    public TruncatedFourierSeries(double[] coefficients, double period) {
        this.coefficients = coefficients;
        this.period = period;
    }

    @Override
    public double value(double x) {
        double result = 0.0;
        int n = coefficients.length;
        if (n % 2 == 0){
            double[] newArray = new double[n + 1];

            // Copy the elements from the original array into the new array
            System.arraycopy(coefficients, 0, newArray, 0, n);

            // Append the new double value at the end of the new array
            newArray[newArray.length - 1] = 0;

            coefficients = newArray;
        }

        // coefficients[0] = A0
        result += coefficients[0];

        // for each set of An, Bn, add the resulting term to the Fourier series
        for (int i = 1; 2*i < coefficients.length; i++) {
            result += coefficients[2*i-1] * Math.cos(2 * Math.PI * i * x / period);
            result += coefficients[2*i] * Math.sin(2 * Math.PI * i * x / period);
        }

        return result;
    }

    /* For debugging */
    public static void main(String[] args) {

        double[] coefficients = {0,1};
        double period = 2 * Math.PI;

        TruncatedFourierSeries fs = new TruncatedFourierSeries(coefficients, period);

        System.out.println("The value is " + fs.value(Math.PI));
    }
}
