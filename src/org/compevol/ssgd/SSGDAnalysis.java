/*
 * SSGDAnalysis.java
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

import dr.evolution.alignment.PatternList;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Arman Bilge
 */
public class SSGDAnalysis implements Bootstrappable {

    private final MaximumLikelihood optimizer;
    private final PairedCompositeLikelihood[] likelihoods;

    public SSGDAnalysis(final MaximumLikelihood optimizer, final PairedCompositeLikelihood... likelihoods) {
        this.optimizer = optimizer;
        this.likelihoods = likelihoods;
    }

    @Override
    public void setPatterns(final PatternList... patterns) {
        // TODO Fix me
//        for (int i = 0; i < likelihoods.length; ++i)
//            likelihoods[i].setPatterns(patterns[i]);
    }

    @Override
    public void run() {
        optimizer.run();
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final MaximumLikelihood optimizer = (MaximumLikelihood) xo.getChild(MaximumLikelihood.class);
            final List<PairedCompositeLikelihood> likelihoods = new ArrayList<PairedCompositeLikelihood>(xo.getChildCount() - 1);
            for (int i = 0; i < xo.getChildCount(); ++i) {
                final Object o = xo.getChild(i);
                if (o instanceof PairedCompositeLikelihood)
                    likelihoods.add((PairedCompositeLikelihood) o);
            }
            return new SSGDAnalysis(optimizer, likelihoods.toArray(new PairedCompositeLikelihood[likelihoods.size()]));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        private final XMLSyntaxRule[] rules = {

        };

        @Override
        public String getParserDescription() {
            return "Performs the SSGD analysis.";
        }

        @Override
        public Class getReturnType() {
            return SSGDAnalysis.class;
        }

        @Override
        public String getParserName() {
            return "ssgd";
        }

    };

}
