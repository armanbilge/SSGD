/*
 * PairedPatterns.java
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

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.math.MathUtils;
import dr.util.Identifiable;

/**
 * @author Arman Bilge <armanbilge@gmail.com>
 */
public class PairedPatterns implements Identifiable {

    private static final long serialVersionUID = 0L;

    private final DataType dataType;
    private final int stateCount;
    private final int transitionCount;
    private final TaxonList taxa;
    private final double[] weights;

    private String id;

    public PairedPatterns(final DataType dataType, final TaxonList taxa) {
        this.dataType = dataType;
        stateCount = dataType.getStateCount();
        transitionCount = stateCount * stateCount;
        this.taxa = taxa;
        weights = new double[taxa.getTaxonCount() * (taxa.getTaxonCount() - 1) / 2 * transitionCount];
    }

    public double getPatternWeight(final Taxon a, final int i, final Taxon b, final int j) {

        final int m = taxa.getTaxonIndex(a);
        final int n = taxa.getTaxonIndex(b);

        if (m == n)
            throw new IllegalArgumentException("The two taxa must be different.");
        else if (m > n)
            return getPatternWeight(b, j, a, i);

        return weights[getIndex(m, n, i, j)];

    }

    public void addPattern(final Taxon a, final int i, final Taxon b, final int j) {
        addPattern(a, i, b, j, 1);
    }

    public void addPattern(final Taxon a, final int i, final Taxon b, final int j, final double w) {

        final int m = taxa.getTaxonIndex(a);
        final int n = taxa.getTaxonIndex(b);

        if (m == n) {
            throw new IllegalArgumentException("The two taxa must be different.");
        } else if (m > n) {
            addPattern(b, j, a, i, w);
            return;
        }

        for (final int x : dataType.getStates(i))
            for (final int y : dataType.getStates(j))
                weights[getIndex(m, n, x, y)] += w;

    }

    private int getIndex(final int m, final int n, final int i, final int j) {
        return (m + n * (n - 1) / 2) * transitionCount + stateCount * i + j;
    }

    public final DataType getDataType() {
        return dataType;
    }

    public final TaxonList getTaxa() {
        return taxa;
    }

    public double[] getFrequencies() {

        final double[] freqs = new double[stateCount];

        for (int i = 0; i < taxa.getTaxonCount(); ++i) {

            for (int j = i+1; j < taxa.getTaxonCount(); ++j) {

                for (int k = 0; k < stateCount; ++k) {

                    for (int l = 0; l < stateCount; ++l) {

                        final double w = weights[getIndex(i, j, k, l)];
                        freqs[k] += w;
                        freqs[l] += w;

                    }

                }
            }

        }

        final double sum = MathUtils.getTotal(freqs);
        for (int i = 0; i < stateCount; ++i)
            freqs[i] /= sum;

        return freqs;
    }

    public double getTotalWeight() {
        return MathUtils.getTotal(weights);
    }

    @Override
    public String getId() {
        return id;
    }

    @Override
    public void setId(final String id) {
        this.id = id;
    }
}
