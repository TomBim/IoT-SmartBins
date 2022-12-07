import pygame

from settings import *
from simulation_page import Simulation_Page
from default_button import Default_Button

DEFAULT_FONT = pygame.font.SysFont('arial', 20)
TITLE_FONT = pygame.font.SysFont('arial', 40)

import bin_lib.some_functions as fcs


class Menu_Page:
    page_name = "menu_page"

    def __init__(self, screen, change_screen, map, everything):
        self.background_color = WHITE
        self.screen = screen

        self.title_surf = TITLE_FONT.render('Simulador SmartBin', False, BLACK)
        self.title_rect = self.title_surf.get_rect(center=(WIDTH//2, 60))
        
        self.see_through = pygame.Surface((WIDTH,HEIGHT)).convert_alpha()
        self.see_through.fill((255, 255, 255, 150))
        self.see_through_rect = self.see_through.get_rect(center=self.screen.get_rect().center)

        self.mapa = map
        self.everything = everything

        self.buttons = [
            Default_Button(screen, WIDTH//2, HEIGHT//2 - 80, 200, 50, DEFAULT_FONT, "Começar simulação",
                           DARK_GRAY, lambda: change_screen(Simulation_Page.page_name)),
            Default_Button(screen, WIDTH//2, HEIGHT//2 +  0, 200, 50, DEFAULT_FONT, "Resetar comércio",
                           DARK_GRAY, lambda: [self.everything.reset_com_points(), fcs.create_rand_com_points(self.mapa, (5,2,4), self.everything)]),
            Default_Button(screen, WIDTH//2, HEIGHT//2 + 80, 200, 50, DEFAULT_FONT, "Sair",
                           DARK_GRAY, lambda: pygame.quit()),
        ]


    def draw_streets(self, screen):
        street_width = 15
        for st in self.mapa.get_streets_list():
            Ax = (st.get_vector()[0].get_pos()[0]) * PROPORTION
            Ay = (st.get_vector()[0].get_pos()[1]) * PROPORTION
            Bx = (st.get_vector()[1].get_pos()[0]) * PROPORTION
            By = (st.get_vector()[1].get_pos()[1]) * PROPORTION
            pygame.draw.line(screen, LIGHT_GRAY, (Ax, Ay), (Bx, By), street_width)


    def draw_intersections(self, screen):
        inter_radius = 7
        for inter in self.mapa.get_intersections_list():
            Ax = (inter.get_pos()[0]) * PROPORTION
            Ay = (inter.get_pos()[1]) * PROPORTION
            pygame.draw.circle(screen, LIGHT_GRAY, (Ax, Ay), inter_radius)


    def draw_comm_point(self, screen):
        comm_radius = 10
        for cp in self.everything._com_points:
            Ax = (cp.get_pos_xy()[0]) * PROPORTION
            Ay = (cp.get_pos_xy()[1]) * PROPORTION
            cp.get_pos_street()
            pygame.draw.circle(screen, DARK_GRAY, (Ax, Ay), comm_radius)
    
    
    def update(self):
        pygame.display.update()
        self.screen.fill(self.background_color)

        self.draw_streets(self.screen)
        self.draw_intersections(self.screen)
        self.draw_comm_point(self.screen)

        self.screen.blit(self.see_through, self.see_through_rect)
        self.screen.blit(self.title_surf, self.title_rect)

        
        if self.buttons:
            for button in self.buttons:
                button.update()
