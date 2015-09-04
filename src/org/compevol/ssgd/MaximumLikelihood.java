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

import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.ml.MLOptimizer;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.ScaleOperator;
import dr.inference.operators.SimpleOperatorSchedule;
import dr.inference.operators.UniformOperator;
import dr.xml.AbstractXMLObjectParser;
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

    private final MLOptimizer optimizer;
    private final Likelihood likelihood;
    private final Parameter variables;

    public MaximumLikelihood(final Likelihood likelihood, final Parameter variables) {
        this.likelihood = likelihood;
        this.variables = variables;
        final SimpleOperatorSchedule schedule = new SimpleOperatorSchedule();
        schedule.addOperator(new ScaleOperator(variables, 100));
        schedule.addOperator(new ScaleOperator(variables, 0.01));
        schedule.addOperator(new UniformOperator(variables, 1.0));
        final Logger[] loggers = {new MCLogger(1)};
        ((MCLogger) loggers[0]).add(likelihood);
        ((MCLogger) loggers[0]).add(variables);
        this.optimizer = new MLOptimizer("ml", 10000, likelihood, schedule, loggers);
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {
        optimizer.run();
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

            final Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
            final Parameter initial = (Parameter) xo.getChild(Parameter.class);

            return new MaximumLikelihood(likelihood, initial);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(Likelihood.class), new ElementRule(Parameter.class)};

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
