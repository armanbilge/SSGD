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

import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class LambertFormatParser extends AbstractXMLObjectParser {

    private static final String FILENAME = "filename";

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final File serializationFile = new File(xo.getId());
        if (serializationFile.exists()) {
            try {
                final ObjectInputStream in = new ObjectInputStream(new FileInputStream(serializationFile));
                try {
                    return in.readObject();
                } catch (final ClassNotFoundException ex) {
                    throw new XMLParseException(ex.getMessage());
                } finally {
                    in.close();
                }
            } catch (final IOException ex) {
                throw new XMLParseException(ex.getMessage());
            }
        }

        final TaxonList taxonList = (TaxonList) xo.getChild(TaxonList.class);

        final Map<String,Taxon> taxa = new HashMap<String,Taxon>(taxonList.getTaxonCount());
        for (int i = 0; i < taxonList.getTaxonCount(); ++i)
            taxa.put(taxonList.getTaxonId(i), taxonList.getTaxon(i));

        final File file = new File(xo.getStringAttribute(FILENAME));
        final PairedPatterns patterns = new PairedPatterns(Nucleotides.INSTANCE, taxonList);

        try {

            final int N;
            {
                int i;
                final BufferedReader br = new BufferedReader(new FileReader(file));
                for (i = 0; br.readLine() != null; ++i);
                br.close();
                N = i;
            }

            for (int i = 0; i < N; ++i) {

                final BufferedReader br = new BufferedReader(new FileReader(file));

                SequenceRecord x = null;
                for (int j = 0; j <= i; ++j) {
                    if (j == i)
                        x = new SequenceRecord(br.readLine());
                    else
                        br.readLine();
                }

                final Taxon a = taxa.get(x.getTaxonName());

                for (int j = i+1; j < N; ++j) {

                    final SequenceRecord y = new SequenceRecord(br.readLine());
                    final Taxon b = taxa.get(y.getTaxonName());

                    patterns.addPattern(a, Nucleotides.A_STATE, b, Nucleotides.A_STATE, x.getACount());
                    patterns.addPattern(a, Nucleotides.C_STATE, b, Nucleotides.C_STATE, x.getCCount());
                    patterns.addPattern(a, Nucleotides.G_STATE, b, Nucleotides.G_STATE, x.getGCount());
                    patterns.addPattern(a, Nucleotides.UT_STATE, b, Nucleotides.UT_STATE, x.getTCount());

                    for (int k = 0; k < x.getSequence().length(); ++k)
                        patterns.addPattern(a, Nucleotides.INSTANCE.getState(x.getSequence().charAt(k)),
                                b, Nucleotides.INSTANCE.getState(y.getSequence().charAt(k)));

                }

                br.close();
            }

        } catch (final IOException ex) {
            throw new XMLParseException(ex.getMessage());
        }

        try {
            final ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream(serializationFile));
            try {
                out.writeObject(patterns);
            } finally {
                out.close();
            }
        } catch (final IOException ex) {
            throw new XMLParseException(ex.getMessage());
        }

        return patterns;
    }

    private static final class SequenceRecord {

        private final String taxon;
        private final long A, T, G, C;
        private final String sequence;

        public SequenceRecord(final String l) {
            final StringTokenizer st = new StringTokenizer(l, ",");
            taxon = st.nextToken();
            A = Long.parseLong(st.nextToken());
            T = Long.parseLong(st.nextToken());
            G = Long.parseLong(st.nextToken());
            C = Long.parseLong(st.nextToken());
            sequence = st.nextToken();
        }

        public String getTaxonName() {
            return taxon;
        }

        public long getACount() {
            return A;
        }

        public long getTCount() {
            return T;
        }

        public long getGCount() {
            return G;
        }

        public long getCCount() {
            return C;
        }

        public String getSequence() {
            return sequence;
        }

    }

    private final XMLSyntaxRule[] rules = {AttributeRule.newStringRule(FILENAME), new ElementRule(TaxonList.class)};

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Converts a file in Lambert format to a PairedPatterns summary.";
    }

    @Override
    public Class getReturnType() {
        return PairedPatterns.class;
    }

    @Override
    public String getParserName() {
        return "lambertFormat";
    }
}
