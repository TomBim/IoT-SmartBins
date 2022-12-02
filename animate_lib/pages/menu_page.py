import pygame

from settings import *
from simulation_page import Simulation_Page
from default_button import Default_Button

DEFAULT_FONT = pygame.font.SysFont('arial', 20)
TITLE_FONT = pygame.font.SysFont('arial', 40)


class Menu_Page:
    page_name = "menu_page"

    def __init__(self, screen, change_screen):
        self.background_color = LIGHT_GRAY
        self.screen = screen
        self.title_surf = TITLE_FONT.render('Simulador SmartBin', False, BLACK)
        self.title_rect = self.title_surf.get_rect(center=(WIDTH//2, 60))

        self.buttons = [
            Default_Button(screen, WIDTH//2, HEIGHT//2 - 40, 200, 50, DEFAULT_FONT, "Começar simulação",
                           DARK_GRAY, lambda: change_screen(Simulation_Page.page_name)),
            Default_Button(screen, WIDTH//2, HEIGHT//2 + 40, 200, 50, DEFAULT_FONT, "Sair",
                           DARK_GRAY, lambda: pygame.quit()),
        ]

    def update(self):
        pygame.display.update()
        self.screen.fill(self.background_color)
        self.screen.blit(self.title_surf, self.title_rect)

        if self.buttons:
            for button in self.buttons:
                button.update()
