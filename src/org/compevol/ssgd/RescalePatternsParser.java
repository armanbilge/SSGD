/*
 * RescalePatternsParser.java
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

import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class RescalePatternsParser extends AbstractXMLObjectParser {

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final List<PairedPatterns> patterns = new ArrayList<PairedPatterns>(xo.getChildCount());
        for (int i = 0; i < xo.getChildCount(); ++i) {
            final Object cxo = xo.getChild(i);
            if (cxo instanceof PairedPatterns)
                patterns.add((PairedPatterns) cxo);
        }

        double sum = 0.0;
        for (final PairedPatterns pattern : patterns)
            sum += pattern.getTotalWeight();

        final double factor = 1.0 / sum;
        Logger.getLogger("org.compevol.ssgd").info("Rescaling all pattern weights by " + factor + ".");
        for (final PairedPatterns pattern : patterns)
            pattern.rescaleWeights(factor);

        return new Object();
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    private final XMLSyntaxRule[] rules = {new ElementRule(PairedPatterns.class, 1, Integer.MAX_VALUE)};

    @Override
    public String getParserDescription() {
        return "rescalePatterns";
    }

    @Override
    public Class getReturnType() {
        return Object.class;
    }

    @Override
    public String getParserName() {
        return "rescalePairedPatterns";
    }
}
