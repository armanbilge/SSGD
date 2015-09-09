/*
 * PairedPatternsSimulator.java
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

import dr.app.seqgen.SeqGen;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import jebl.evolution.alignments.Alignment;
import jebl.evolution.sequences.Sequence;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class PairedPatternsSimulator {

    private final CoalescentSimulator coalescentSimulator = new CoalescentSimulator();
    private final SeqGen sequenceSimulator;

    private final TaxonList taxa;
    private final DemographicModel demographicModel;
    private final int locusCount;

    public PairedPatternsSimulator(final TaxonList taxa, final DemographicModel demographicModel, final SiteModel siteModel, final int locusLength, final int locusCount) {
        this.taxa = taxa;
        this.demographicModel = demographicModel;
        sequenceSimulator = new SeqGen(locusLength, 1.0, siteModel.getFrequencyModel(), siteModel.getSubstitutionModel(), siteModel, 0.0);
        this.locusCount = locusCount;
    }

    private Tree simulateTree() {
        return coalescentSimulator.simulateTree(taxa, demographicModel);
    }

    private Alignment simulateSite() {
        return sequenceSimulator.simulate(simulateTree());
    }

    public PairedPatterns simulatePatterns() {

        final PairedPatterns patterns = new PairedPatterns(Nucleotides.INSTANCE, taxa);

        for (int i = 0; i < locusCount; ++i) {

            final Alignment alignment = simulateSite();

            for (int j = 0; j < taxa.getTaxonCount(); ++j) {

                final Taxon a = taxa.getTaxon(j);
                final Sequence x = alignment.getSequence(jebl.evolution.taxa.Taxon.getTaxon(a.getId()));

                for (int k = j+1; k < taxa.getTaxonCount(); ++k) {

                    final Taxon b = taxa.getTaxon(k);
                    final Sequence y = alignment.getSequence(jebl.evolution.taxa.Taxon.getTaxon(b.getId()));

                    for (int l = 0; l < alignment.getSiteCount(); ++l)
                        patterns.addPattern(a, x.getState(l).getIndex(), b, y.getState(l).getIndex());

                }

            }

        }

        return patterns;

    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String LENGTH = "length";
        private static final String LOCI = "loci";

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            return new PairedPatternsSimulator(
                    (TaxonList) xo.getChild(TaxonList.class),
                    (DemographicModel) xo.getChild(DemographicModel.class),
                    (SiteModel) xo.getChild(SiteModel.class),
                    xo.getIntegerAttribute(LENGTH),
                    xo.getIntegerAttribute(LOCI)).simulatePatterns();
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }
        private final XMLSyntaxRule[] rules = {
                new ElementRule(TaxonList.class),
                new ElementRule(DemographicModel.class),
                new ElementRule(SiteModel.class),
                AttributeRule.newIntegerRule(LENGTH),
                AttributeRule.newIntegerRule(LOCI)
        };

        @Override
        public String getParserDescription() {
            return "Simulates a PairedPatterns object.";
        }

        @Override
        public Class getReturnType() {
            return PairedPatterns.class;
        }

        @Override
        public String getParserName() {
            return "pairedPatternsSimulator";
        }

    };

}
