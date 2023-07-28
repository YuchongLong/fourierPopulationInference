package beast.base.evolution.tree.coalescent;


import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.integration.IterativeLegendreGaussIntegrator;
import org.apache.commons.math3.analysis.integration.LegendreGaussIntegrator;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import static org.apache.commons.math.analysis.solvers.BrentSolver.DEFAULT_ABSOLUTE_ACCURACY;
import static org.apache.commons.math3.analysis.integration.BaseAbstractUnivariateIntegrator.*;

@Description("coalescent intervals for a Fourier population")

public class FourierPopulation extends PopulationFunction.Abstract {
    final public Input<Function> coeffParameter = new Input<>("coefficients",
            "Fourier coefficients", Validate.REQUIRED);
    final public Input<RealParameter> periodParameter = new Input<>("periodicity",
            "periodicity of the periodic functions", Validate.REQUIRED);


    // Implementation of abstract methods

    @Override
	public List<String> getParameterIds() {
    	List<String> ids = new ArrayList<>();
    	if (coeffParameter.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)coeffParameter.get()).getID());
        if (periodParameter.get() instanceof BEASTInterface)
            ids.add(((BEASTInterface)periodParameter.get()).getID());
        return ids;
    }


    @Override
	public double getPopSize(double t) {
        // logN = A0 + sum(An*cos(n*pi*t)+Bn*sin(n*pi*t)), so N = exp(A0 + sum(An*cos(n*t)+Bn*sin(n*t)))
        UnivariateFunction logFourierPopSize = new TruncatedFourierSeries(coeffParameter.get().getDoubleValues(), periodParameter.get().getValue());
        return Math.exp(logFourierPopSize.value(t));
    }


    @Override
	public double getIntensity(double t) {

        // Create an instance of the LegendreGaussIntegrator
        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(10,DEFAULT_RELATIVE_ACCURACY,
                DEFAULT_ABSOLUTE_ACCURACY,DEFAULT_MIN_ITERATIONS_COUNT,DEFAULT_MAX_ITERATIONS_COUNT);

        UnivariateFunction fourierPopSize = new TruncatedFourierSeries(coeffParameter.get().getDoubleValues(), periodParameter.get().getValue());

        UnivariateFunction expNegFourierPopSize = x -> Math.exp(-fourierPopSize.value(x));

        // Perform the integration
        if (t == 0) {
            return 0;
        }else {
            return integrator.integrate(Integer.MAX_VALUE, expNegFourierPopSize, 0, t); // will throw an error if t=0
        }
    }

    @Override
	public double getInverseIntensity(double x) {
        //throw new UnsupportedOperationException("getInverseIntensity not implemented");
        return 2; // or any number larger than 1
    }


    /* For debugging */
    public static void main(String[] args) {

        /*
        <distribution spec="Coalescent" tree="@tree">
          <populationModel spec="FourierPopulation">
            <coefficients spec="RealParameter" value="0.1 14.3 12"/>
          </populationModel>
        </distribution>
         */
        FourierPopulation fp = new FourierPopulation();
        fp.initByName("coefficients", new RealParameter("1 1 1"), "periodicity", new RealParameter("Math.PI"));
//        fp.initByName("coefficients", new RealParameter("0.243 0.262 1.817E32"));

        // System.out.println("Pop size at time t=3.0: " + fp.getPopSize(3.0));

        // Write population function variation through time to a file

        try (PrintStream ps = new PrintStream("file_path")) {
            ps.println("t N intensity");

            double T = 10;
            double steps = 101;
            double dt = T/(steps - 1);

            for (int i=0; i<steps; i++) {
                double t = i*dt;
                ps.println(t + " " + fp.getPopSize(t) + " " + fp.getIntensity(t));
            }

        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        }

//        // Test IterativeLegendreGaussIntegrator //
//        IterativeLegendreGaussIntegrator integrator = new IterativeLegendreGaussIntegrator(10,DEFAULT_RELATIVE_ACCURACY,
//                DEFAULT_ABSOLUTE_ACCURACY,DEFAULT_MIN_ITERATIONS_COUNT,DEFAULT_MAX_ITERATIONS_COUNT);
//
//        UnivariateFunction fn = x -> x+Math.pow(x,2);
//
//        System.out.println("Integration result = " + integrator.integrate(Integer.MAX_VALUE, fn,0, 1));
    }


}
