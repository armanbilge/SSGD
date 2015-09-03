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
import lbfgsb.Bound;
import lbfgsb.DifferentiableFunction;
import lbfgsb.LBFGSBException;
import lbfgsb.Minimizer;
import lbfgsb.Result;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class MaximumLikelihood implements Spawnable {

    private final Minimizer optimizer;
    private final DifferentiableFunction likelihood;
    private final Parameter variables;
    private final double[] initial;

    public MaximumLikelihood(final DifferentiableFunction likelihood, final Parameter variables) {
        this.likelihood = likelihood;
        this.variables = variables;
        initial = new double[variables.getDimension()];
        Arrays.fill(initial, 1);
//        this.initial = variables.getParameterValues();
        this.optimizer = new Minimizer();
        final List<Bound> bounds = new ArrayList<Bound>(variables.getDimension());
        final Bounds<Double> myBounds = variables.getBounds();
        for (int i = 0; i < myBounds.getBoundsDimension(); ++i)
            bounds.add(new Bound(myBounds.getLowerLimit(i) / variables.getValue(i), myBounds.getUpperLimit(i) / variables.getValue(i)));
        optimizer.setBounds(bounds);
        optimizer.setCorrectionsNo(100);
        optimizer.getStopConditions().setFunctionReductionFactor(1.0);
        optimizer.getStopConditions().setFunctionReductionFactorInactive();
        optimizer.getStopConditions().setMaxGradientNormInactive();
        optimizer.getStopConditions().setMaxIterationsInactive();
        optimizer.setDebugLevel(2);
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {
        try {
            final Result result = optimizer.run(likelihood, initial);
            System.out.println(result);
        } catch (final LBFGSBException ex) {
            throw new RuntimeException(ex);
        }
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

            final DifferentiableFunction likelihood = (DifferentiableFunction) xo.getChild(DifferentiableFunction.class);
            final Parameter initial = (Parameter) xo.getChild(Parameter.class);

            return new MaximumLikelihood(likelihood, initial);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(DifferentiableFunction.class), new ElementRule(Parameter.class)};

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
