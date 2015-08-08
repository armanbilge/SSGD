/*
 * LambertFormatParser.java
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
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxa;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Arrays;
import java.util.InputMismatchException;
import java.util.Scanner;
import java.util.logging.Logger;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class LambertFormatParser extends AbstractXMLObjectParser {

    private static final Logger LOGGER = Logger.getLogger("org.compevol.ssgd");
    private static final String FILENAME = "filename";

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final Taxa taxa = (Taxa) xo.getChild(Taxa.class);

        final int[] constantA = new int[taxa.getTaxonCount()];
        Arrays.fill(constantA, Nucleotides.INSTANCE.getState('A'));
        final int[] constantC = new int[taxa.getTaxonCount()];
        Arrays.fill(constantC, Nucleotides.INSTANCE.getState('C'));
        final int[] constantG = new int[taxa.getTaxonCount()];
        Arrays.fill(constantG, Nucleotides.INSTANCE.getState('G'));
        final int[] constantT = new int[taxa.getTaxonCount()];
        Arrays.fill(constantT, Nucleotides.INSTANCE.getState('T'));

        final File file = new File(xo.getStringAttribute(FILENAME));
        final Patterns patterns = new Patterns(Nucleotides.INSTANCE, taxa);

        try {

            final Scanner sc = new Scanner(file);

            while (sc.hasNextLine()) {

                final SequenceRecord[] records = new SequenceRecord[taxa.getTaxonCount()];
                for (int i = 0; i < records.length; ++i) {
                    records[i] = new SequenceRecord(sc.nextLine());
                }

                LOGGER.info("Processing gene " + records[0].getGeneName() + ".");

                patterns.addPattern(constantA, records[0].getACount());
                patterns.addPattern(constantC, records[0].getCCount());
                patterns.addPattern(constantG, records[0].getGCount());
                patterns.addPattern(constantT, records[0].getTCount());

                for (int i = 0; i < records[0].getSequence().length(); ++i) {
                    final int[] pattern = new int[taxa.getTaxonCount()];
                    for (final SequenceRecord record : records) {
                        final int j = taxa.getTaxonIndex(record.getTaxonName());
                        if (j == -1)
                            throw new XMLParseException("No taxon with id=" + record.getTaxonName());
                        pattern[j] = Nucleotides.INSTANCE.getState(record.getSequence().charAt(i));
                    }
                    patterns.addPattern(pattern, 1.0);
                }

            }

            sc.close();

        } catch (final FileNotFoundException ex) {
            throw new XMLParseException(ex.getMessage());
        }

        return patterns;
    }

    private static int nextInt(final Scanner sc) {
        try {
            return sc.nextInt();
        } catch (final InputMismatchException ime) {
            return 0;
        }
    }

    private static final class SequenceRecord {

        private final String gene;
        private final String taxon;
        private final int A, T, G, C;
        private final String sequence;

        public SequenceRecord(final String l) {
            final Scanner sc = new Scanner(l);
            sc.useDelimiter(",");
            gene = sc.next();
            taxon = sc.next();
            A = nextInt(sc);
            T = nextInt(sc);
            G = nextInt(sc);
            C = nextInt(sc);
            sequence = sc.next();
        }

        public String getGeneName() {
            return gene;
        }

        public String getTaxonName() {
            return taxon;
        }

        public int getACount() {
            return A;
        }

        public int getTCount() {
            return T;
        }

        public int getGCount() {
            return G;
        }

        public int getCCount() {
            return C;
        }

        public String getSequence() {
            return sequence;
        }

    }

    private final XMLSyntaxRule[] rules = {AttributeRule.newStringRule(FILENAME), new ElementRule(Taxa.class)};

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Converts a file in Lambert format to a BEAST PatternList.";
    }

    @Override
    public Class getReturnType() {
        return PatternList.class;
    }

    @Override
    public String getParserName() {
        return "lambertFormat";
    }
}
