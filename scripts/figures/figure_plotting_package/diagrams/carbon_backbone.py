from .config import Vector, ParameterName, ColorConfig, ZOrderConfig
from .config import CompositeFigure, PathStep, PathOperation, PathShape, Rectangle, Circle, cos_sin


class CarbonBackboneConfig(object):
    common_z_order = ZOrderConfig.default_patch_z_order
    z_order_increment = ZOrderConfig.z_order_increment
    ear2_z_order = common_z_order + z_order_increment

    carbon_radius = 0.04
    unlabeled_color = ColorConfig.light_blue
    labeled_color = ColorConfig.slightly_light_orange

    single_carbon_backbone_config = {
        ParameterName.edge_width: None,
        ParameterName.radius: carbon_radius,
        ParameterName.z_order: common_z_order
    }


class SingleCarbonBackbone(CompositeFigure):
    def __init__(
            self, center: Vector, carbon_num: int, radius: float, angle=0,
            scale=1, bottom_left_offset=None, base_z_order=0, z_order_increment=1, face_color=None, **kwargs):
        short_edge = radius * 2
        long_edge = radius * 2 * carbon_num
        assert -90 <= angle <= 90
        if -45 <= angle <= 45:
            width = long_edge
            height = short_edge
            bound_angle = angle
        else:
            if angle > 45:
                bound_angle = angle - 90
            else:
                bound_angle = angle + 90
            width = short_edge
            height = long_edge
        bound_rectangle = Rectangle(center, width, height, bound_angle)
        # circle_center_dict = {
        #     f'circle_{index}': center - radius * (2 * index - carbon_num + 1) * cos_sin(angle)
        #     for index in range(carbon_num)}
        circle_name_list = [f'circle_{index}' for index in range(carbon_num)]
        circle_center_list = [
            center - radius * (2 * index - carbon_num + 1) * cos_sin(angle) for index in range(carbon_num)]
        if isinstance(face_color, (tuple, list)):
            assert len(face_color) == carbon_num
            circle_color_list = face_color
        else:
            circle_color_list = [face_color] * carbon_num
        # circle_obj_dict = {
        #     name: Circle(**{
        #         ParameterName.center: circle_center,
        #         ParameterName.radius: radius,
        #         ParameterName.name: name,
        #         ParameterName.z_order: ZOrderConfig.default_patch_z_order,
        #         **kwargs
        #     })
        #     for name, circle_center in circle_center_dict.items()
        # }
        circle_obj_dict = {
            name: Circle(**{
                ParameterName.center: circle_center,
                ParameterName.radius: radius,
                ParameterName.name: name,
                ParameterName.z_order: ZOrderConfig.default_patch_z_order,
                ParameterName.face_color: circle_color,
                **kwargs
            })
            for name, circle_center, circle_color in zip(circle_name_list, circle_center_list, circle_color_list)
        }
        single_carbon_backbone_dict = {
            'circles': circle_obj_dict,
        }
        super().__init__(
            single_carbon_backbone_dict, bottom_left=bound_rectangle.bottom_left, size=bound_rectangle.size,
            scale=scale, bottom_left_offset=bottom_left_offset, base_z_order=base_z_order,
            z_order_increment=z_order_increment)


class CarbonBackbone(CompositeFigure):
    total_width = 1
    height_to_width_ratio = 1

    def __init__(
            self, carbon_num: int, labeled=True, single_carbon_backbone=False, color_dict=None,
            scale=1, bottom_left_offset=None, base_z_order=0, z_order_increment=1, **kwargs):
        total_width = self.total_width
        total_height = total_width * self.height_to_width_ratio
        if color_dict is None:
            color_dict = (
                CarbonBackboneConfig.labeled_color if labeled else CarbonBackboneConfig.unlabeled_color)
        if single_carbon_backbone:
            carbon_backbone_center_angle_dict_list = [
                {
                    ParameterName.center: Vector(0.5, 0.5),
                },
            ]
        else:
            carbon_backbone_center_angle_dict_list = [
                {
                    ParameterName.center: Vector(0.65, 0.15),
                    ParameterName.angle: 30,
                },
                {
                    ParameterName.center: Vector(0.4, 0.8),
                    ParameterName.angle: 10,
                },
                {
                    ParameterName.center: Vector(0.2, 0.5),
                    ParameterName.angle: 23,
                },
                {
                    ParameterName.center: Vector(0.31, 0.23),
                    ParameterName.angle: -10,
                },
                {
                    ParameterName.center: Vector(0.55, 0.61),
                    ParameterName.angle: 41,
                },
            ]
        carbon_backbone_common_dict = {
            **CarbonBackboneConfig.single_carbon_backbone_config,
            ParameterName.carbon_num: carbon_num,
            ParameterName.face_color: color_dict,
            **kwargs
        }

        background_box = Rectangle(**{
            ParameterName.center: Vector(0.5, 0.5),
            ParameterName.width: total_width,
            ParameterName.height: total_height,
            ParameterName.face_color: ColorConfig.light_gray,
            ParameterName.z_order: 0
        })
        single_carbon_backbone_obj_dict = {
            f'single_backbone_{index}':
                SingleCarbonBackbone(**{
                    **carbon_backbone_common_dict,
                    **single_backbone_center_angle_dict,
                    ParameterName.name: f'single_backbone_{index}',
                }) for index, single_backbone_center_angle_dict in enumerate(carbon_backbone_center_angle_dict_list)
        }
        carbon_backbone_dict = {
            # ParameterName.background: {'background': background_box},
            'single_carbon_backbone': single_carbon_backbone_obj_dict,
        }
        super().__init__(
            carbon_backbone_dict, bottom_left=Vector(0, 0), size=Vector(total_width, total_height),
            scale=scale, bottom_left_offset=bottom_left_offset, base_z_order=base_z_order,
            z_order_increment=z_order_increment)


class ScatteredLabeledCarbonBackbone(CompositeFigure):
    total_width = 0.8
    height_to_width_ratio = 0.55

    @staticmethod
    def labeled_color_converter(labeled_indicator_list, color_dict):
        # return [
        #     CarbonBackboneConfig.labeled_color if labeled_indicator else CarbonBackboneConfig.unlabeled_color
        #     for labeled_indicator in labeled_indicator_list]
        return [color_dict[labeled_indicator] for labeled_indicator in labeled_indicator_list]

    def __init__(
            self, carbon_num: int=0, single_carbon_backbone=False, selected_carbon_list=(),
            color_dict=None, scale=1, bottom_left_offset=None, base_z_order=0, z_order_increment=1, **kwargs):
        total_width = self.total_width
        total_height = total_width * self.height_to_width_ratio
        if color_dict is None:
            color_dict = {1: CarbonBackboneConfig.labeled_color, 0: CarbonBackboneConfig.unlabeled_color}
        if single_carbon_backbone:
            carbon_backbone_center_angle_dict_list = [
                {
                    ParameterName.center: Vector(total_width, total_height) / 2,
                    ParameterName.carbon_num: carbon_num,
                    ParameterName.face_color: self.labeled_color_converter(selected_carbon_list, color_dict),
                },
            ]
        else:
            carbon_backbone_center_angle_dict_list = [
                {
                    ParameterName.center: Vector(0.4, 0.35),
                    ParameterName.carbon_num: 6,
                    ParameterName.face_color: self.labeled_color_converter([1, 0, 1, 0, 0, 1], color_dict),
                },
                {
                    ParameterName.center: Vector(0.15, 0.18),
                    ParameterName.carbon_num: 3,
                    ParameterName.face_color: self.labeled_color_converter([0, 1, 1], color_dict),
                },
                {
                    ParameterName.center: Vector(0.25, 0.06),
                    ParameterName.carbon_num: 3,
                    ParameterName.face_color: self.labeled_color_converter([1, 1, 0], color_dict),
                },
                {
                    ParameterName.center: Vector(0.62, 0.11),
                    ParameterName.carbon_num: 3,
                    ParameterName.face_color: self.labeled_color_converter([1, 1, 1], color_dict),
                },
                {
                    ParameterName.center: Vector(0.5, 0.23),
                    ParameterName.carbon_num: 4,
                    ParameterName.face_color: self.labeled_color_converter([1, 1, 0, 0], color_dict),
                },
            ]
        carbon_backbone_common_dict = {
            **CarbonBackboneConfig.single_carbon_backbone_config,
            **kwargs
        }

        single_carbon_backbone_obj_dict = {
            f'single_backbone_{index}':
                SingleCarbonBackbone(**{
                    **carbon_backbone_common_dict,
                    **single_backbone_center_angle_dict,
                    ParameterName.name: f'single_backbone_{index}',
                }) for index, single_backbone_center_angle_dict in enumerate(carbon_backbone_center_angle_dict_list)
        }
        carbon_backbone_dict = {
            'single_carbon_backbone': single_carbon_backbone_obj_dict,
        }
        super().__init__(
            carbon_backbone_dict, bottom_left=Vector(0, 0), size=Vector(total_width, total_height),
            scale=scale, bottom_left_offset=bottom_left_offset, base_z_order=base_z_order,
            z_order_increment=z_order_increment, background=False,)
