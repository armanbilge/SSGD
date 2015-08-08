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

import dr.math.DifferentialEvolution;
import dr.math.MachineAccuracy;
import dr.math.MultivariateFunction;
import dr.math.MultivariateMinimum;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.Spawnable;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class MaximumLikelihood implements Spawnable {

    private final MultivariateMinimum optimizer;
    private final MultivariateFunction likelihood;
    private final double[] initial;

    public MaximumLikelihood(final MultivariateFunction likelihood, final double[] initial) {
        this.likelihood = likelihood;
        this.initial = initial;
        optimizer = new DifferentialEvolution(likelihood.getNumArguments());
    }

    public MaximumLikelihood(final MultivariateFunction likelihood, final double[] initial, final int populationSize) {
        this.likelihood = likelihood;
        this.initial = initial;
        optimizer = new DifferentialEvolution(likelihood.getNumArguments(), populationSize);
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {
        optimizer.optimize(new NegativeMultivariateFunction(likelihood), initial, MachineAccuracy.EPSILON, MachineAccuracy.EPSILON);
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String INITIAL = "initial";
        private static final String POPULATION_SIZE = "populationSize";

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

            final MultivariateFunction likelihood = (MultivariateFunction) xo.getChild(MultivariateFunction.class);
            final double[] initial = xo.getDoubleArrayAttribute(INITIAL);

            if (xo.hasAttribute(POPULATION_SIZE))
                return new MaximumLikelihood(likelihood, initial, xo.getIntegerAttribute(POPULATION_SIZE));
            else
                return new MaximumLikelihood(likelihood, initial);
            
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(MultivariateFunction.class),
                AttributeRule.newDoubleArrayRule(INITIAL), AttributeRule.newIntegerRule(POPULATION_SIZE, true)};

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
