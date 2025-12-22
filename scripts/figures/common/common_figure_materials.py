from .config import plt, mcolors, ColorConfig


class ColorMapName(object):
    diet_color_map = 'diet_color_map'
    tissue_color_map = 'tissue_color_map'
    state_color_map = 'state_color_map'
    blue_white_orange_cmap = 'blue_white_orange_cmap'


set_2_color = plt.get_cmap('Set2').colors
pastel_1_color = plt.get_cmap('Pastel1').colors


class CommonColorDict(object):
    pass


class CommonFigureString(object):
    pass
