import time

import pygame

from text_box import Text_box


class Stats:
    """Track statistics for Gravity Simulator."""

    STATSBOARD_FONT_SIZE = 20
    STATSBOARD_SIZE_X = 260
    STATSBOARD_SIZE_Y = 23
    FIXED_STEP_SIZE_INTEGRATORS_COLOR = (120, 233, 250)
    ADAPTIVE_STEP_SIZE_INTEGRATORS_COLOR = (173, 255, 47)

    def __init__(self, grav_sim) -> None:
        self.simulation_time = 0
        self.fps = grav_sim.clock.get_fps()
        self.total_energy = 0
        self.settings = grav_sim.settings
        self.is_paused = False
        self.is_holding_rclick = False
        self._create_statsboard(grav_sim)
        self._statsboard_init_print_msg()

    def update(self, grav_sim) -> None:
        self.fps = grav_sim.clock.get_fps()
        self.objects_count = len(grav_sim.grav_objs)

        if grav_sim.menu.main_menu_active == True:
            self.start_time = time.time()
        if self.is_paused == False:
            self.run_time = time.time() - self.start_time

        # Count time for the create new star function
        if self.is_holding_rclick == True:
            self.holding_rclick_time = time.time() - self.holding_rclick_start_time

    def reset(self, grav_sim) -> None:
        self.start_time = time.time()
        self.simulation_time = 0
        self.total_energy = 0
        grav_sim.simulator.is_initialize = True
        grav_sim.simulator.is_initialize_integrator = (
            grav_sim.simulator.current_integrator
        )
        grav_sim.camera._pos[0] = 0
        grav_sim.camera._pos[1] = 0

    def start_pause(self) -> None:
        self.paused_start_time = time.time()
        self.is_paused = True

    def end_pause(self) -> None:
        self.start_time -= self.paused_start_time - time.time()
        self.is_paused = False

    def start_holding_rclick(self) -> None:
        self.holding_rclick_start_time = time.time()
        self.is_holding_rclick = True

    def end_holding_rclick(self) -> None:
        self.is_holding_rclick = False

    def print_msg(self) -> None:
        self.fps_board.print_msg(f"FPS = {self.fps:2.1f}")
        self.obj_count_board.print_msg(f"Object = {self.objects_count}")
        self.simulation_time_board.print_msg(
            f"Simulation Time = {self.simulation_time / 365.242189:.1e} years"
        )
        self.run_time_board.print_msg(f"Run time = {self.run_time:.0f} seconds")
        self.total_energy_board.print_msg(f"Total Energy = {self.total_energy:.3e}")

        self.star_img_scale_board.print_msg(
            f"Star Image Scale = {self.settings.star_img_scale}"
        )
        self.planet_img_scale_board.print_msg(
            f"Planet Image Scale = {self.settings.planet_img_scale}"
        )
        self.distance_scale_board.print_msg(
            f"Distance Scale = {self.settings.distance_scale}"
        )
        self.new_star_mass_scale_board.print_msg(
            f"New star mass scale = {self.settings.new_star_mass_scale:g}x"
        )
        self.new_star_speed_scale_board.print_msg(
            f"New star speed scale = {self.settings.new_star_speed_scale:d}x"
        )
        self.dt_board.print_msg(f"dt = {self.settings.dt:g} days / frame")
        self.time_speed_board.print_msg(f"Time Speed = {self.settings.time_speed}x")
        self.max_iteration_board.print_msg(
            f"Max iterations / frame = {self.settings.max_iteration}"
        )
        self.min_iteration_board.print_msg(
            f"Min iterations / frame = {self.settings.min_iteration}"
        )
        self.tolerance_board.print_msg(f"Tolerance = {self.settings.tolerance:g}")

    def draw(self, grav_sim) -> None:
        self.print_msg()
        self.fps_board.draw()
        self.obj_count_board.draw()
        self.simulation_time_board.draw()
        self.run_time_board.draw()
        self.total_energy_board.draw()

        self.parameters_board.draw()
        self.star_img_scale_board.draw()
        self.planet_img_scale_board.draw()
        self.distance_scale_board.draw()
        self.new_star_mass_scale_board.draw()
        self.new_star_speed_scale_board.draw()
        self.dt_board.draw()
        self.time_speed_board.draw()
        self.max_iteration_board.draw()
        self.min_iteration_board.draw()
        self.tolerance_board.draw()

        self.integrators_board.draw()
        self.fixed_step_size_board.draw()
        self.euler_board.draw()
        self.euler_cromer_board.draw()
        self.rk4_board.draw()
        self.leapfrog_board.draw()

        self.adaptive_step_size_board.draw()
        self.rkf45_board.draw()
        self.dopri_board.draw()
        self.dverk_board.draw()
        self.rkf78_board.draw()
        self.ias15_board.draw()

        # Visual indicator for currently changing parameter
        match self.settings.current_changing_parameter:
            case "star_img_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.star_img_scale_board.rect.centery + 5),
                    4,
                )
            case "planet_img_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.planet_img_scale_board.rect.centery + 5),
                    4,
                )
            case "distance_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.distance_scale_board.rect.centery + 5),
                    4,
                )
            case "new_star_mass_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.new_star_mass_scale_board.rect.centery + 5),
                    4,
                )
            case "new_star_speed_scale":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.new_star_speed_scale_board.rect.centery + 5),
                    4,
                )
            case "dt":
                pygame.draw.circle(
                    grav_sim.screen, "yellow", (290, self.dt_board.rect.centery + 5), 4
                )
            case "time_speed":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.time_speed_board.rect.centery + 5),
                    4,
                )
            case "max_iteration":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.max_iteration_board.rect.centery + 5),
                    4,
                )
            case "min_iteration":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.min_iteration_board.rect.centery + 5),
                    4,
                )
            case "tolerance":
                pygame.draw.circle(
                    grav_sim.screen,
                    "yellow",
                    (290, self.tolerance_board.rect.centery + 5),
                    4,
                )

        # Visual indicator for currently selected integrator
        match grav_sim.simulator.current_integrator:
            case "euler":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.euler_board.rect.centery + 5),
                    4,
                )
            case "euler_cromer":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.euler_cromer_board.rect.centery + 5),
                    4,
                )
            case "rk4":
                pygame.draw.circle(
                    grav_sim.screen, "green", (290, self.rk4_board.rect.centery + 5), 4
                )
            case "leapfrog":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.leapfrog_board.rect.centery + 5),
                    4,
                )
            case "rkf45":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.rkf45_board.rect.centery + 5),
                    4,
                )
            case "dopri":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.dopri_board.rect.centery + 5),
                    4,
                )
            case "dverk":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.dverk_board.rect.centery + 5),
                    4,
                )
            case "rkf78":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.rkf78_board.rect.centery + 5),
                    4,
                )
            case "ias15":
                pygame.draw.circle(
                    grav_sim.screen,
                    "green",
                    (290, self.ias15_board.rect.centery + 5),
                    4,
                )

    def check_button(self, grav_sim, mouse_pos) -> None:
        """Check if there is any click on the buttons"""
        if self.settings.is_hide_gui == False:
            if self.star_img_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_star_img_scale = True
            if self.planet_img_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_planet_img_scale = True
            if self.distance_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_distance_scale = True
            if self.new_star_mass_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_new_star_mass_scale = True
            if self.new_star_speed_scale_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_new_star_speed_scale = True
            if self.dt_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_dt = True
            if self.time_speed_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_time_speed = True
            if self.max_iteration_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_max_iteration = True
            if self.min_iteration_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_min_iteration = True
            if self.tolerance_board.rect.collidepoint(mouse_pos):
                self.settings.set_all_parameters_changing_false()
                self.settings.is_changing_tolerance = True

            if self.euler_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_euler = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "euler"
            if self.euler_cromer_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_euler_cromer = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "euler_cromer"
            if self.rk4_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rk4 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rk4"
            if self.leapfrog_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_leapfrog = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "leapfrog"

            if self.rkf45_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rkf45 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rkf45"
            if self.dopri_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_dopri = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "dopri"
            if self.dverk_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_dverk = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "dverk"
            if self.rkf78_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_rkf78 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "rkf78"
            if self.ias15_board.rect.collidepoint(mouse_pos):
                grav_sim.simulator.set_all_integrators_false()
                grav_sim.simulator.is_ias15 = True
                grav_sim.simulator.is_initialize = True
                grav_sim.simulator.is_initialize_integrator = "ias15"

    def _statsboard_init_print_msg(self) -> None:
        self.parameters_board.print_msg("Parameters: (Click below to select)")
        self.integrators_board.print_msg("Integrators: (Click below to select)")
        self.fixed_step_size_board.print_msg("(Fixed Step Size)")
        self.euler_board.print_msg("Euler")
        self.euler_cromer_board.print_msg("Euler-Cromer")
        self.rk4_board.print_msg("4th order Runge-Kutta")
        self.leapfrog_board.print_msg("Leapfrog (Verlet)")
        self.adaptive_step_size_board.print_msg("(Adaptive Step Size)")
        self.rkf45_board.print_msg("Runge-Kutta-Fehleberg 4(5)")
        self.dopri_board.print_msg("Dormand-Prince 5(4)")
        self.dverk_board.print_msg("Verner's method 6(5) DVERK")
        self.rkf78_board.print_msg("Runge-Kutta-Fehlberg 7(8)")
        self.ias15_board.print_msg("IAS15")

    @classmethod
    def _create_statsboard(self, grav_sim) -> None:
        self.fps_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 0),
        )
        self.obj_count_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 23),
        )

        self.simulation_time_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 46),
        )
        self.run_time_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 69),
        )
        self.total_energy_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 92),
        )

        self.parameters_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 138),
        )
        self.star_img_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 161),
        )
        self.planet_img_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 184),
        )
        self.distance_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 207),
        )
        self.new_star_mass_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 230),
        )
        self.new_star_speed_scale_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 253),
        )
        self.dt_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 276),
            text_color=self.FIXED_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.time_speed_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 299),
            text_color=self.FIXED_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.max_iteration_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 322),
            text_color=self.ADAPTIVE_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.min_iteration_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 345),
            text_color=self.ADAPTIVE_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.tolerance_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 368),
            text_color=self.ADAPTIVE_STEP_SIZE_INTEGRATORS_COLOR,
        )

        self.integrators_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 414),
        )
        self.fixed_step_size_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 437),
            text_color=self.FIXED_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.euler_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 460),
        )
        self.euler_cromer_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 483),
        )
        self.rk4_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 506),
        )
        self.leapfrog_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 529),
        )
        self.adaptive_step_size_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 575),
            text_color=self.ADAPTIVE_STEP_SIZE_INTEGRATORS_COLOR,
        )
        self.rkf45_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 598),
        )
        self.dopri_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 621),
        )
        self.dverk_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 644),
        )
        self.rkf78_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 667),
        )
        self.ias15_board = Text_box(
            grav_sim,
            self.STATSBOARD_FONT_SIZE,
            size_x=self.STATSBOARD_SIZE_X,
            size_y=self.STATSBOARD_SIZE_Y,
            font="Manrope",
            text_box_left_top=(10, 713),
        )
