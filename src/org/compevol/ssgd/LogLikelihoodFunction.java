/*
 * LogLikelihoodFunction.java
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

import com.cureos.numerics.Calcfc;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class LogLikelihoodFunction implements Calcfc {

    private final Likelihood function;
    private final Parameter variables;
    private final double[] scale;

    public LogLikelihoodFunction(final Likelihood function, final Parameter variables) {
        this.function = function;
        this.variables = variables;
        scale = variables.getAttributeValue();
    }

    @Override
    public double Compute(final int n, final int m, final double[] x, final double[] con) {
        for (int i = 0; i < variables.getDimension(); ++i)
            variables.setParameterValue(i, x[i] * scale[i]);
        for (int i = 0, j = variables.getDimension(); i < variables.getDimension(); ++i, ++j) {
            con[i] = x[i] - variables.getBounds().getLowerLimit(i) / scale[i];
            con[j] = variables.getBounds().getUpperLimit(i) / scale[i] - x[i];
        }
        return -function.getLogLikelihood();
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            return new LogLikelihoodFunction((Likelihood) xo.getChild(Likelihood.class), (Parameter) xo.getChild(Parameter.class));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        final XMLSyntaxRule[] rules = {new ElementRule(Likelihood.class), new ElementRule(Parameter.class)};

        @Override
        public String getParserDescription() {
            return "A multivariate function view for a likelihood.";
        }

        @Override
        public Class getReturnType() {
            return LogLikelihoodFunction.class;
        }

        @Override
        public String getParserName() {
            return "logLikelihoodFunction";
        }
    };
}
