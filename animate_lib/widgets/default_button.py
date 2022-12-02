import pygame

class Default_Button:
    def __init__(self, screen, x, y, w, h, font, text, color, function):
        self.screen = screen
        self.function = function 
        self.position = (x,y)
        
        # text_font = pygame.font.Font(PIXELED_FONT, size)
        self.text_surf = font.render(text, True, (0, 0, 0))
        self.text_rect = self.text_surf.get_rect(center = self.position)
        
        self.button_bg_image = pygame.Surface((w, h))
        self.button_bg_image.fill(color)
        self.button_bg_rect = self.button_bg_image.get_rect(center = self.position)

        self.was_pressed = False


    def move_button(self):
        self.button_bg_rect.center = (self.position[0], self.position[1] - 3)
        self.text_rect.center = (self.position[0], self.position[1] - 3)


    def move_button_back(self):
        self.button_bg_rect.center = self.position
        self.text_rect.center = self.position


    def update(self):
        self.screen.blit(self.button_bg_image, self.button_bg_rect)
        self.screen.blit(self.text_surf, self.text_rect)
            
        if (not pygame.mouse.get_pressed()[0]) and self.button_bg_rect.collidepoint(pygame.mouse.get_pos()):
            if self.was_pressed: # was released
                self.function()
        
        if self.button_bg_rect.collidepoint(pygame.mouse.get_pos()): # is hovering
            self.move_button()
            if pygame.mouse.get_pressed()[0]: # was pressed
                self.was_pressed = True
                self.move_button_back()
            else:
                self.was_pressed = False
        else:
            self.move_button_back()
            self.was_pressed = False