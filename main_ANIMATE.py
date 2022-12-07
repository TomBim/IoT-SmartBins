import pygame, sys, os
pygame.font.init()

sys.path.append('./bin_lib')
sys.path.append('./animate_lib')
sys.path.append('./animate_lib/pages')
sys.path.append('./animate_lib/widgets')

from settings import *
from menu_page import Menu_Page
from simulation_page import Simulation_Page

import bin_lib.map as map
import bin_lib.entities as entities
import bin_lib.some_functions as fcs

FILE_INTERSECTIONS = "bin_lib/intersections.txt"
FILE_STRETS = "bin_lib/streets.txt"
FILE_COMMERTIAL_POINTS = ""


class Simulation:
    def __init__(self):
        pygame.init()
        
        pygame.display.set_caption('Simulador Smart Bin')
        self.screen = pygame.display.set_mode((WIDTH, HEIGHT))
        pygame.font.init()

        
        self.clock = pygame.time.Clock()

        self.mapa = map.read_map(FILE_INTERSECTIONS, FILE_STRETS, FILE_COMMERTIAL_POINTS)
        self.everything = entities.Everything(self.mapa)
        fcs.create_rand_com_points(self.mapa, (5,2,4), self.everything)

        self.screen_name = Menu_Page.page_name
        self.page = Menu_Page(self.screen,self.change_screen, self.mapa, self.everything)

        

    def run(self):
        while(True):
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    pygame.quit()
                    sys.exit()
            
            self.page.update()
            self.clock.tick(FPS)


    def change_screen(self, screen_name):
        self.screen_name = screen_name
        self.page = {
            Menu_Page.page_name : Menu_Page(self.screen, self.change_screen, self.mapa, self.everything),
            Simulation_Page.page_name : Simulation_Page(self.screen, self.change_screen, self.mapa, self.everything),
        }.get(screen_name, Menu_Page.page_name)


    def go_to_Menu_page(self):
        self.change_screen(Menu_Page.page_name)


if __name__ == '__main__':
    simulation = Simulation()
    simulation.run()