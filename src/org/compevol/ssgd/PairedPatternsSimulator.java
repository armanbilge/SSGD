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
import dr.evomodel.coalescent.PiecewisePopulationModel;
import dr.evomodel.sitemodel.SiteModel;
import jebl.evolution.alignments.Alignment;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class PairedPatternsSimulator {

    private final CoalescentSimulator coalescentSimulator = new CoalescentSimulator();
    private final SeqGen sequenceSimulator;

    private final TaxonList taxa;
    private final PiecewisePopulationModel populationModel;
    private final int length;

    public PairedPatternsSimulator(final TaxonList taxa, final PiecewisePopulationModel populationModel, final SiteModel siteModel, final int length) {
        this.taxa = taxa;
        this.populationModel = populationModel;
        sequenceSimulator = new SeqGen(1, 1.0, siteModel.getFrequencyModel(), siteModel.getSubstitutionModel(), siteModel, 0.0);
        this.length = length;
    }

    private Tree simulateTree() {
        return coalescentSimulator.simulateTree(taxa, populationModel);
    }

    private Alignment simulateSite() {
        return sequenceSimulator.simulate(simulateTree());
    }

    public PairedPatterns simulatePatterns() {

        final PairedPatterns patterns = new PairedPatterns(Nucleotides.INSTANCE, taxa);

        for (int i = 0; i < length; ++i) {

            final Alignment site = simulateSite();

            for (int j = 0; j < taxa.getTaxonCount(); ++j) {

                for (int k = j+1; k < taxa.getTaxonCount(); ++k) {

                    final Taxon a = taxa.getTaxon(j);
                    final Taxon b = taxa.getTaxon(k);
                    final int x = site.getSequence(jebl.evolution.taxa.Taxon.getTaxon(a.getId())).getState(0).getIndex();
                    final int y = site.getSequence(jebl.evolution.taxa.Taxon.getTaxon(b.getId())).getState(0).getIndex();
                    patterns.addPattern(a, x, b, y);

                }

            }

        }

        return patterns;

    }

}
