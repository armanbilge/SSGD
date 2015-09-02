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
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.BOBYQAOptimizer;

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
        this.initial = variables.getParameterValues();
        optimizer = new BOBYQAOptimizer(variables.getDimension() * 3 / 2);
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {
        final double[] lower = new double[variables.getDimension()];
        final double[] upper = new double[variables.getDimension()];
        final Bounds<Double> bounds = variables.getBounds();
        for (int i = 0; i < bounds.getBoundsDimension(); ++i) {
            lower[i] = bounds.getLowerLimit(i);
            upper[i] = bounds.getUpperLimit(i);
        }
        final PointValuePair maximum = optimizer.optimize(new ObjectiveFunction(likelihood), new InitialGuess(initial), new SimpleBounds(lower, upper), MaxEval.unlimited(), GoalType.MAXIMIZE);
        System.out.println(Arrays.toString(maximum.getPoint()));
        System.out.println(maximum.getValue());
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
