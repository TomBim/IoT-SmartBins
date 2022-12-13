import pygame
import random

from settings import *
from default_button import Default_Button

DEFAULT_FONT = pygame.font.SysFont('arial', 20)

# import bin_lib.map as map
# import bin_lib.entities as entities
import bin_lib.some_functions as fcs

TIME_STEP = 1
# bins/streets are emptied/cleaned once each:
FRAME_TO_EMPTY_BINS = 60*30
FRAME_TO_CLEAN_STREETS = 60*20


def draw_health_bar(surf, pos, size, borderC, backC, healthC, progress):
    pygame.draw.rect(surf, backC, (*pos, *size))
    pygame.draw.rect(surf, borderC, (*pos, *size), 1)
    innerPos  = (pos[0]+1, pos[1]+1)
    innerSize = ((size[0]-2) * progress, size[1]-2)
    rect = (round(innerPos[0]), round(innerPos[1]), round(innerSize[0]), round(innerSize[1]))
    pygame.draw.rect(surf, healthC, rect)


class Simulation_Page:
    page_name = "simulation_page"

    def __init__(self, screen, change_screen, map, everything):
        self.background_color = WHITE
        self.screen = screen

        self.mapa = map
        self.everything = everything

        self.everything.reset()
        fcs.create_rand_bins(self.mapa, self.everything)
        self.frame_count_sweep_streets = 0
        self.frame_count_empty_bins = 0

        self.TRASH_BIN = pygame.image.load('animate_lib/assets/trash_bin.png').convert_alpha()
        self.TRASH_BIN = pygame.transform.scale(self.TRASH_BIN, (20,24))

        self.view_path = True

        self.buttons = [
            Default_Button(screen, 50, 50, 40, 40, DEFAULT_FONT, "x", DARK_GRAY, lambda: change_screen("menu_page")),
            Default_Button(screen, 220, 50, 250, 40, DEFAULT_FONT, "Mostrar/Esconder Rotas", DARK_GRAY, lambda: self.toggle_view_path()),
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
                BLACK, WHITE, RED, percentage)

    def draw_trash(self, screen):
        trash_radius = 2
        for t in self.everything._pos_trash_floor:
            Ax = (t.get_pos_xy()[0]) * PROPORTION - 5 + 10*random.random()
            Ay = (t.get_pos_xy()[1]) * PROPORTION - 5 + 10*random.random()
            pygame.draw.circle(screen, BLACK, (Ax, Ay), trash_radius)

    def draw_people(self, screen):
        people_radius = 4
        for p in self.everything._ppl:
            Ax = (p.get_pos_xy()[0]) * PROPORTION
            Ay = (p.get_pos_xy()[1]) * PROPORTION
            # If impatient, paint RED
            if ((p._max_time_carrying_trash - p._time_carrying_trash) > 5) or not p.has_trash():
                color = BLUE
            else: color = RED
            # Draw person 
            pygame.draw.circle(screen, color, (Ax, Ay), people_radius)
            # Draw person fov when with trash
            if p._has_trash:
                pygame.draw.circle(screen, color, (Ax, Ay), p.get_fov() * PROPORTION, 2)
            if self.view_path:
                # Draw person destination
                destinX = p.get_destination().get_pos_xy()[0] * PROPORTION
                destinY = p.get_destination().get_pos_xy()[1] * PROPORTION
                pygame.draw.line(screen, GREEN, (Ax, Ay), (destinX, destinY), 1)
                # Draw person next intersection
                if len(p.get_path()) > 0:
                    pathX = p.get_path()[0].get_pos()[0] * PROPORTION
                    pathY = p.get_path()[0].get_pos()[1] * PROPORTION
                    pygame.draw.line(screen, BLACK, (Ax, Ay), (pathX, pathY), 1)

    
    def draw_trash_counter(self, screen):
        value = self.everything._trash_in_the_streets
        value = f'{value:.2f} L'
        screen.blit(DEFAULT_FONT.render(str(value), True, BLACK), (WIDTH-100, 50))


    def toggle_view_path(self):
        self.view_path = not self.view_path
            
    def update(self):
        pygame.display.update()
        self.screen.fill(self.background_color)

        self.draw_streets(self.screen)
        self.draw_intersections(self.screen)
        self.draw_comm_point(self.screen)
        
        self.draw_bins(self.screen)
        self.draw_trash(self.screen)
        self.draw_people(self.screen)
        self.draw_trash_counter(self.screen)

        self.everything.update_people(TIME_STEP, 0)
        fcs.create_rand_ppl(self.mapa, self.everything, TIME_STEP)
        if self.frame_count_sweep_streets // FRAME_TO_CLEAN_STREETS:
            self.everything.sweep_streets()
            self.frame_count_sweep_streets = 0
        
        if (self.frame_count_empty_bins % FRAME_TO_EMPTY_BINS) == 0 and self.frame_count_empty_bins > 0:
            self.everything.empty_bins(self.frame_count_empty_bins)
        
        self.frame_count_sweep_streets += 1
        self.frame_count_empty_bins += 1

        if self.buttons:
            for button in self.buttons:
                button.update()

