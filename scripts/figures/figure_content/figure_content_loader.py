from ..common.config import enum
from ..figure_plotting_package.main import draw_figure


class BaseFigureName(enum.Enum):
    def __str__(self):
        return self.value

    def __eq__(self, other):
        return self.value.__eq__(other)

    @staticmethod
    def return_all():
        pass


def test_figure_content_loader(figure_direct, figure_name):
    import importlib
    current_file_name = f'figure_{figure_name}'
    current_direct = __name__[:__name__.rindex('.')]
    current_file_path = f'{current_direct}.{figure_direct}.{current_file_name}'
    imported_lib = importlib.import_module(current_file_path)
    Figure = imported_lib.Figure
    return Figure


class NormalFigureName(BaseFigureName):
    figure_4 = '1'
    figure_s4 = 's4'
    figure_5 = '5'
    figure_s5 = 's5'
    test_figures = 'test_figures'

    @staticmethod
    def return_test_figures():
        from .figures.test_figures import load_all_figure_obj
        return load_all_figure_obj()


FigureName = NormalFigureName


def normal_figure_content_loader(figure_name):
    if figure_name == NormalFigureName.figure_4:
        from .figures.figure_4 import Figure
    elif figure_name == NormalFigureName.figure_s4:
        from .figures.figure_s4 import Figure
    elif figure_name == NormalFigureName.figure_5:
        from .figures.figure_5 import Figure
    elif figure_name == NormalFigureName.figure_s5:
        from .figures.figure_s5 import Figure
    else:
        Figure = test_figure_content_loader('figures', figure_name)
    return Figure()


def figure_plotting_main(figure_name, output_svg=False):
    if figure_name == NormalFigureName.test_figures:
        figure_obj_list = NormalFigureName.return_test_figures()
        for figure_obj in figure_obj_list:
            draw_figure(figure_obj, output_svg)
    else:
        figure_obj = normal_figure_content_loader(figure_name)
        draw_figure(figure_obj, output_svg)
