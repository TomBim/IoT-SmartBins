import pygame, sys, os
pygame.font.init()

sys.path.append('./bin_lib')
sys.path.append('./animate_lib')
sys.path.append('./animate_lib/pages')
sys.path.append('./animate_lib/widgets')

from settings import *
from menu_page import Menu_Page
from simulation_page import Simulation_Page


class Simulation:
    def __init__(self):
        pygame.init()
        
        pygame.display.set_caption('Simulador Smart Bin')
        self.screen = pygame.display.set_mode((WIDTH, HEIGHT))
        
        self.clock = pygame.time.Clock()
        
        self.screen_name = Menu_Page.page_name
        self.page = Menu_Page(self.screen,self.change_screen)


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
            Menu_Page.page_name : Menu_Page(self.screen, self.change_screen),
            Simulation_Page.page_name : Simulation_Page(self.screen, self.change_screen)
        }.get(screen_name, Menu_Page.page_name)


    def go_to_Menu_page(self):
        self.change_screen(Menu_Page.page_name)


if __name__ == '__main__':
    simulation = Simulation()
    simulation.run()