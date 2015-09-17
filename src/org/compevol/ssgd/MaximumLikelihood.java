/*
 * Estimate.java
 *
 * SSGD: Serially-Sampled Genome Demographics
 *
 * Copyright (c) 2015 Arman Bilge <armanbilge@gmail.com>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package org.compevol.ssgd;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MachineAccuracy;
import dr.math.MathUtils;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.Spawnable;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.SimpleBounds;
import org.apache.commons.math3.optim.SimplePointChecker;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.CMAESOptimizer;
import org.apache.commons.math3.random.RandomGenerator;

import java.util.Arrays;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class MaximumLikelihood implements Spawnable {

    private final MultivariateOptimizer optimizer;
    private final MultivariateFunction likelihood;
    private final Parameter variables;
    private final double[] initial;

    public MaximumLikelihood(final MultivariateFunction likelihood, final Parameter variables) {
        this.likelihood = likelihood;
        this.variables = variables;
        initial = new double[variables.getDimension()];
        Arrays.fill(initial, 1);
        optimizer = new CMAESOptimizer(Integer.MAX_VALUE, 0.0, true, 0, 8096, new RandomGenerator() {

            @Override
            public void setSeed(int i) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void setSeed(int[] ints) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void setSeed(long l) {
                throw new UnsupportedOperationException();
            }

            @Override
            public void nextBytes(byte[] bytes) {
                MathUtils.nextBytes(bytes);
            }

            @Override
            public int nextInt() {
                return MathUtils.nextInt();
            }

            @Override
            public int nextInt(int i) {
                return MathUtils.nextInt(i);
            }

            @Override
            public long nextLong() {
                return MathUtils.nextLong();
            }

            @Override
            public boolean nextBoolean() {
                return MathUtils.nextBoolean();
            }

            @Override
            public float nextFloat() {
                return MathUtils.nextFloat();
            }

            @Override
            public double nextDouble() {
                return MathUtils.nextDouble();
            }

            @Override
            public double nextGaussian() {
                return MathUtils.nextGaussian();
            }
        }, true, new SimplePointChecker<PointValuePair>(MachineAccuracy.SQRT_EPSILON, MachineAccuracy.EPSILON));
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {
        final Bounds<Double> bounds = variables.getBounds();
        final double[] lower = new double[bounds.getBoundsDimension()];
        final double[] upper = new double[bounds.getBoundsDimension()];
        final double[] sigma = new double[variables.getDimension()];
        Arrays.fill(sigma, 1.0);
        for (int i = 0; i < variables.getDimension(); ++i) {
            lower[i] = bounds.getLowerLimit(i) / variables.getParameterValue(i);
            upper[i] = bounds.getUpperLimit(i) / variables.getParameterValue(i);
        }
        final PointValuePair result = optimizer.optimize(
                new CMAESOptimizer.PopulationSize(4 + 3 * (int) Math.log(variables.getDimension())),
                new CMAESOptimizer.Sigma(sigma),
                GoalType.MAXIMIZE,
                new ObjectiveFunction(likelihood),
                new InitialGuess(initial),
                new SimpleBounds(lower, upper),
                new MaxEval(Integer.MAX_VALUE)
        );
        System.out.println(variables);
        System.out.println(result.getValue());
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

            final MultivariateFunction likelihood = (MultivariateFunction) xo.getChild(MultivariateFunction.class);
            final Parameter initial = (Parameter) xo.getChild(Parameter.class);

            return new MaximumLikelihood(likelihood, initial);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(MultivariateFunction.class), new ElementRule(Parameter.class)};

        @Override
        public String getParserDescription() {
            return "Maximizes the arguments for the given likelihood function.";
        }

        @Override
        public Class getReturnType() {
            return MaximumLikelihood.class;
        }

        @Override
        public String getParserName() {
            return "maximumLikelihood";
        }
    };

}
