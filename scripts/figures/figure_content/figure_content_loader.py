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
    figure_1 = '1'
    fangchao_fly_mid_figures = 'fangchao_fly_mid'

    @staticmethod
    def return_fangchao_fly_mid_figures():
        from .figures.figure_fangchao_fly_mid import load_all_figure_obj
        return load_all_figure_obj()


FigureName = NormalFigureName


def normal_figure_content_loader(figure_name):
    if figure_name == NormalFigureName.figure_1:
        from .figures.figure_1 import Figure
    else:
        Figure = test_figure_content_loader('figures', figure_name)
    return Figure()


def figure_plotting_main(figure_name, output_svg=False):
    if figure_name == NormalFigureName.fangchao_fly_mid_figures:
        figure_obj_list = NormalFigureName.return_fangchao_fly_mid_figures()
        for figure_obj in figure_obj_list:
            draw_figure(figure_obj, output_svg)
    else:
        figure_obj = normal_figure_content_loader(figure_name)
        draw_figure(figure_obj, output_svg)
