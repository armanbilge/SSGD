/*
 * SSGD.java
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

import dr.app.plugin.Plugin;
import dr.xml.XMLObjectParser;

import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class SSGD implements Plugin {

    private final Set<XMLObjectParser> parsers;

    {
        final Set<XMLObjectParser> parsers = new HashSet<XMLObjectParser>();
        parsers.add(SSGDAnalysis.PARSER);
        parsers.add(PairedCompositeLikelihood.PARSER);
        parsers.add(HKYSkylineIntegrator.PARSER);
        parsers.add(TaxonSpecificSequenceErrorModel.PARSER);
        parsers.add(new LambertFormatParser());
        parsers.add(PairedPatternsSimulator.PARSER);
        parsers.add(new PairedPatternsFrequenciesParser());
        parsers.add(MaximumLikelihood.PARSER);
        parsers.add(LogLikelihoodFunction.PARSER);
        parsers.add(Bootstrapper.PARSER);
        this.parsers = Collections.unmodifiableSet(parsers);
    }

    @Override
    public Set<XMLObjectParser> getParsers() {
        return parsers;
    }

}
