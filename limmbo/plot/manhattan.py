

from numpy import arange, asarray, cumsum, flipud, log10, logical_and


def plot_manhattan(df,
                   alpha=None,
                   null_style=dict(alpha=0.1, color='Orange'),
                   alt_style=dict(alpha=1.0, color='DarkBlue'),
                   ax=None,
                   annotate=False):
    r"""Produce a manhattan plot.

    Arguments:
        df : :class:`pandas.DataFrame`
            A Pandas DataFrame containing columns pv for p-values, pos for
            base-pair positions, and chrom for chromossome names..
        alpha : float
            Threshold for significance. Defaults to 0.01 significance level
            (bonferroni-adjusted).
        ax : :class:`matplotlib.axes.AxesSubplot`:
            The target handle for this figure. If None, the current axes is
            set.

    Returns:
        (:class:`matplotlib.axes.AxesSubplot`)
            Axes object.

    Examples:

    .. plot::
        from numpy.random import RandomState
        from numpy import arange, ones, kron
        import pandas as pd
        from limix.plot import plot_manhattan
        from matplotlib import pyplot as plt
        random = RandomState(1)
        pv = random.rand(5000)
        pv[1200:1250] = random.rand(50)**4
        chrom  = kron(arange(1,6), ones(1000))
        pos = kron(ones(5), arange(1,1001))
        data = dict(pv=pv, chrom=chrom, pos=pos)
        plot_manhattan(pd.DataFrame(data=data))
        plt.tight_layout()
        plt.show()
    """

    import matplotlib.pyplot as plt

    ax = plt.gca() if ax is None else ax
    df['chrom'] = df['chrom'].astype(int)

    if 'pos' not in df:
        df['pos'] = arange(df.shape[0])
    else:
        df['pos'] = df['pos'].astype(int)

    df = df.sort_values(['chrom', 'pos'])
    if 'label' not in df:
        chrom = df['chrom'].astype(int).astype(str)
        pos = df['pos'].astype(int).astype(str)
        df['label'] = (
            'chrom' + chrom + '_pos' + pos)

    df = _abs_pos(df)

    if alpha is None:
        alpha = 0.01 / df.shape[0]

    ytop = -1.2 * log10(min(df['pv'].min(), alpha))

    _plot_chrom_strips(ax, df, ytop)
    _plot_points(ax, df, alpha, null_style, alt_style, annotate=annotate)
    _set_frame(ax, df, ytop)

    ax.set_ylabel('-log$_{10}$pv')
    ax.set_xlabel('chromosome')

    _set_ticks(ax, _chrom_bounds(df), df['chrom'].unique())

    return ax


def _set_frame(ax, df, ytop):
    ax.set_ylim(0, ytop)
    ax.set_xlim(0, df['abs_pos'].max())

    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)


def _plot_points(ax, df, alpha, null_style, alt_style, annotate):
    null_df = df.loc[df['pv'] >= alpha, :]
    alt_df = df.loc[df['pv'] < alpha, :]

    ax.plot(null_df['abs_pos'], -log10(null_df['pv']), '.', ms=5, **null_style)
    ax.plot(alt_df['abs_pos'], -log10(alt_df['pv']), '.', ms=5, **alt_style)
    if annotate:
        for i in range(alt_df.shape[0]):
            x = alt_df['abs_pos'].values[i]
            y = -log10(alt_df['pv'].values[i])
            _annotate(ax, x, y, alt_df['label'].values[i])


def _plot_chrom_strips(ax, df, ytop):
    uchroms = df['chrom'].unique()
    chrom_bounds = _chrom_bounds(df)
    for i in range(0, len(uchroms), 2):
        ax.fill_between(
            x=chrom_bounds,
            y1=0,
            y2=ytop,
            where=logical_and(chrom_bounds >= chrom_bounds[i],
                              chrom_bounds <= chrom_bounds[i + 1]),
            facecolor='LightGray',
            linewidth=0,
            alpha=0.5)


def _set_ticks(ax, chrom_bounds, uchroms):
    n = len(uchroms)
    xticks = asarray([chrom_bounds[i:i + 2].mean() for i in range(n)])
    ax.set_xticks(xticks)
    ax.tick_params(axis='x', which='both', labelsize=6)
    ax.set_xticklabels(uchroms)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')


def _abs_pos(df):
    uchroms = df['chrom'].unique()
    chrom_ends = [int(df['pos'][df['chrom'] == c].max()) for c in uchroms]

    offset = flipud(cumsum(chrom_ends)[:-1])

    df['abs_pos'] = df['pos'].copy()

    uchroms = list(reversed(uchroms))
    for i in range(len(offset)):
        ix = df['chrom'] == uchroms[i]
        df.loc[ix, 'abs_pos'] = df.loc[ix, 'abs_pos'] + offset[i]

    return df


def _chrom_bounds(df):
    uchroms = df['chrom'].unique()
    min_pos = [df['abs_pos'][df['chrom'] == c].min() for c in uchroms]
    max_pos = [df['abs_pos'][df['chrom'] == c].max() for c in uchroms]
    v = []
    for i in range(len(uchroms) - 1):
        v.append((max_pos[i] + min_pos[i + 1]) / 2)
    return asarray([min_pos[0]] + v + [max_pos[len(v)]])


def _annotate(ax, x, y, text):
    ax.annotate(
        text,
        xy=(x, y),
        xytext=(-18, 18),
        textcoords='offset points',
        fontsize=6,
        ha='center',
        va='bottom',
        bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.3),
        arrowprops=dict(
           arrowstyle='->', connectionstyle='arc3,rad=0.5', color='red'))
