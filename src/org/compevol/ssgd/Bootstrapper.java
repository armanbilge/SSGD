/*
 * Bootstrapper.java
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
import dr.inference.loggers.Logger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.Spawnable;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import jdk.nashorn.internal.runtime.linker.Bootstrap;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Arman Bilge
 */
public class Bootstrapper implements Spawnable {

    private final Bootstrappable analysis;
    private final Logger logger;
    private final int replicates;
    private final PatternList[] patterns;

    public Bootstrapper(final Bootstrappable analysis, final Logger logger, final int replicates, final PatternList... patterns) {
        this.analysis = analysis;
        this.logger = logger;
        this.replicates = replicates;
        this.patterns = patterns;
    }

    @Override
    public boolean getSpawnable() {
        return true;
    }

    @Override
    public void run() {

        logger.startLogging();

        analysis.setPatterns(patterns);
        analysis.run();
        logger.log(0);

        for (int i = 1; i <= replicates; ++i) {

            final PatternList[] bootstrappedPatterns = new PatternList[patterns.length];
            for (int j = 0; j < patterns.length; ++j)
                bootstrappedPatterns[j] = new BootstrappedPatterns(patterns[j]);

            analysis.setPatterns(bootstrappedPatterns);
            analysis.run();

            logger.log(i);

        }

    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String REPLICATES = "replicates";

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

            final Bootstrappable analysis = (Bootstrappable) xo.getChild(Bootstrappable.class);
            final Logger logger = (Logger) xo.getChild(Logger.class);
            final int replicates = xo.getIntegerAttribute(REPLICATES);
            final List<PatternList> patterns = new ArrayList<PatternList>();
            for (int i = 0; i < xo.getChildCount(); ++i) {
                final Object o = xo.getChild(i);
                if (o instanceof PatternList)
                    patterns.add((PatternList) o);
            }

            return new Bootstrapper(analysis, logger, replicates, patterns.toArray(new PatternList[patterns.size()]));
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        private final XMLSyntaxRule[] rules = {new ElementRule(Bootstrappable.class),
                new ElementRule(Logger.class),
                AttributeRule.newIntegerRule(REPLICATES),
                new ElementRule(PatternList.class, 1, Integer.MAX_VALUE)};

        @Override
        public String getParserDescription() {
            return "Runs and bootstraps an analysis.";
        }

        @Override
        public Class getReturnType() {
            return Bootstrap.class;
        }

        @Override
        public String getParserName() {
            return "bootstrapper";
        }

    };

}
