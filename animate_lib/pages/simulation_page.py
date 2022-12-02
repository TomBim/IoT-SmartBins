import pygame

from settings import *
from default_button import Default_Button

DEFAULT_FONT = pygame.font.SysFont('arial', 20)

import bin_lib.map as map
import bin_lib.entities as entities
import bin_lib.some_functions as fcs


FILE_INTERSECTIONS = "bin_lib/intersections.txt"
FILE_STRETS = "bin_lib/streets.txt"
FILE_COMMERTIAL_POINTS = ""

TIME_OF_SIMULATION = 3600*365*5
TIME_STEP = 1

PROPORTION = WIDTH/540


def draw_health_bar(surf, pos, size, borderC, backC, healthC, progress):
    pygame.draw.rect(surf, backC, (*pos, *size))
    pygame.draw.rect(surf, borderC, (*pos, *size), 1)
    innerPos  = (pos[0]+1, pos[1]+1)
    innerSize = ((size[0]-2) * progress, size[1]-2)
    rect = (round(innerPos[0]), round(innerPos[1]), round(innerSize[0]), round(innerSize[1]))
    pygame.draw.rect(surf, healthC, rect)


class Simulation_Page:
    page_name = "simulation_page"

    def __init__(self, screen, change_screen):
        self.background_color = WHITE
        self.screen = screen

        self.mapa = map.read_map(FILE_INTERSECTIONS, FILE_STRETS, FILE_COMMERTIAL_POINTS)
        self.everything = entities.Everything(self.mapa)
        fcs.create_rand_com_points(self.mapa, (5,2,4), self.everything)
        fcs.create_rand_bins(self.mapa, self.everything)

        self.TRASH_BIN = pygame.image.load('animate_lib/assets/trash_bin.png').convert_alpha()
        self.TRASH_BIN = pygame.transform.scale(self.TRASH_BIN, (20,24))

        self.buttons = [
            Default_Button(screen, 50, 50 , 40, 40, DEFAULT_FONT, "x", (150, 150, 150), lambda: change_screen("menu_page")),
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
        inter_radius = 5 
        for inter in self.mapa.get_intersections_list():
            Ax = (inter.get_pos()[0]) * PROPORTION
            Ay = (inter.get_pos()[1]) * PROPORTION
            pygame.draw.circle(screen, LIGHT_GRAY, (Ax, Ay), inter_radius)


    def draw_comm_point(self, screen):
        comm_radius = 10
        for cp in self.everything._com_points:
            Ax = (cp.get_pos_xy()[0]) * PROPORTION
            Ay = (cp.get_pos_xy()[1]) * PROPORTION
            pygame.draw.circle(screen, DARK_GRAY, (Ax, Ay), comm_radius)


    def draw_bins(self, screen):
        bin_radius = 10
        for b in self.everything._bins:
            Ax = (b.get_pos_xy()[0]) * PROPORTION
            Ay = (b.get_pos_xy()[1]) * PROPORTION
            rect = self.TRASH_BIN.get_rect()
            rect.center = (Ax, Ay)
            screen.blit(self.TRASH_BIN, rect)
            self.draw_health(rect, screen, b.get_percentage())

        
    def draw_health(self, rect, surf, percentage):
        health_rect = pygame.Rect(0, 0, self.TRASH_BIN.get_width(), 7)
        health_rect.midbottom = rect.centerx, rect.top
        max_health = 10
        draw_health_bar(surf, health_rect.topleft, health_rect.size, 
                (0, 0, 0), (255, 0, 0), (0, 255, 0), percentage)


    def draw_people(self, screen):
        people_radius = 2
        for p in self.everything._ppl:
            Ax = (p.get_pos_xy()[0]) * PROPORTION
            Ay = (p.get_pos_xy()[1]) * PROPORTION
            pygame.draw.circle(screen, BLUE, (Ax, Ay), people_radius)
            pygame.draw.circle(screen, BLUE, (Ax, Ay), p.get_fov() * PROPORTION, 2)
            
            
    def update(self):
        pygame.display.update()
        self.screen.fill(self.background_color)

        self.draw_streets(self.screen)
        self.draw_intersections(self.screen)
        
        self.draw_comm_point(self.screen)
        self.draw_bins(self.screen)
        self.draw_people(self.screen)

        self.everything.update_people(TIME_STEP)
        fcs.create_rand_ppl(self.mapa, self.everything, TIME_STEP)

        if self.buttons:
            for button in self.buttons:
                button.update()
