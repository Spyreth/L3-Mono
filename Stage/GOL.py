import numpy as np
import pygame
import time
import sys

"""
CONWAY'S GAME OF LIFE : PYGAME IMPLEMENTATION

To run :  python GOL.py (width of a cell in pixels) (number of cells on axis x) (number of cells on axis y) (fps_max)
ex : python GOL.py 10 80 60 60

CONTROLS :
- Spacebar : pause or restart the game  /!\ THE GAME IS NOT RUNNING WHEN STARTED
- Rightclick : Revert the state of a cell (alive -> dead , dead -> alive)


Implemented :
- Control over nb of cells and cells width
- Possibility to stop the game
- Possibility to draw alive and dead cells

To do :
- Color control
- Check for MOUSE_BUTTON_UP to make the drawing easier
- Grid control
- Able/disable the periodicity of the grid
- Make the program run faster or slower without modifying the ability to draw cells fast
"""

colorB = (10,10,10)  # dead cells color (black)
colorW = (255,255,255)  # live cells color (white)
colorGrid = (40,40,40)  # grid color FOR IMPLEMENTATION



def update_cells(cells):
    """Given an array of cells, generates the next one using the rules of the Game Of Life
    Args:
        cells (array): 2D numpy array containing the state of each cell at t0
    Returns:
        updated_cells (array): 2D numpy array containing the state of each cell at t1
    """

    cells_updated = np.zeros((cells.shape[0],cells.shape[1]))  #Initialization with only zeros of the new array

    #Loop over every cell in the array
    for i in range(cells.shape[0]):
        for j in range(cells.shape[1]):

            neighbours = np.sum(cells[i-1:i+2 , j-1:j+2]) - cells[i,j]  #Calcul of the number of neighbors for each cell
            
            if cells[i,j] == 1 :
                if 2 <= neighbours <= 3:    #if alive and 2/3 neighbours -> alive
                    cells_updated[i,j] = 1
            else :
                if neighbours == 3:   #if dead and 3 neighbours -> alive
                    cells_updated[i,j] = 1 
            #for every other case the cell is/becomes dead, we keep the value 0
                
    return cells_updated



def update_screen(screen, updated_cells, size):
    """Creates the image in pygame for the given cell repartition

    Args:
        screen (pygame.screen): The screen object from pygame
        updated_cells (array): array of cells to generate the image
        size (int): size in pixels of one cell

    Returns:
        void
    """

    for i in range(updated_cells.shape[0]):

        for j in range(updated_cells.shape[1]):

            if updated_cells[i,j] == 0:
                color = color = colorB
            else:
                color = colorW
            
            pygame.draw.rect(screen, color, (i*size,j*size,size,size))

    return 0



def update(screen,cells,size):
    """Updates the screen from its previous state. Combination of the previous 2 funcs to only loop one time
    Args:
        cells (array): 2D numpy array containing the state of each cell at t0
        screen (pygame.screen): The screen object from pygame
        size (int): size in pixels of one cell
    Returns:
        cells_updated (array): 2D numpy array containing the state of each cell at t1
    """

    cells_updated = np.zeros((cells.shape[0],cells.shape[1]))  #Initialization with only zeros of the new array

    #Loop over every cell in the array
    for i in range(cells.shape[0]):
        for j in range(cells.shape[1]):

            neighbours = np.sum(cells[i-1:i+2 , j-1:j+2]) - cells[i,j]  #Calcul of the number of neighbors
            color = colorB
            
            if cells[i,j] == 1 :
                if 2 <= neighbours <= 3:    #if alive and 2/3 neighbours -> alive
                    cells_updated[i,j] = 1
                    color = colorW
            else :
                if neighbours == 3:   #if dead and 3 neighbours -> alive
                    cells_updated[i,j] = 1
                    color = colorW
            #for every other case the cell is/becomes dead, we keep the value 0 and color black
            
            pygame.draw.rect(screen, color, (i*size,j*size,size,size)) #pygame drawing of the cell
                
    return cells_updated
    
    
    

def main(cell_width, dim1, dim2, fps_max):
    """Main function of the program, makes the game run by updating the cell and the screen and checking for pygame envents.

    Args:
        cell_width (int): size of one cell in pixel
        dim1 (int): nb of cells in the x axis
        dim2 (int): nb of cells in the y axis
    """

    #initialization of pygame
    pygame.init()
    pygame.display.set_caption("conway's game of life")
    clock = pygame.time.Clock()

    screendim1 = dim1*cell_width
    screendim2 = dim2*cell_width

    screen = pygame.display.set_mode((screendim1,screendim2))

    #Cells + screen initialization
    cells0 = np.zeros((dim1,dim2))
    #cells2 = (np.random.rand(dim1,dim2)<0.4).astype(int)    #random start, needs imlementationi
    screen.fill(colorB)
    cells = update(screen, cells0, cell_width)
    pygame.display.flip()

    running = False  #The cells are not evolving at first

    #game loop
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
                
            elif event.type == pygame.KEYDOWN:
                if event.key == pygame.K_SPACE:
                    running = not running
                    cells = update(screen,cells,cell_width)
                    pygame.display.flip()
            if pygame.mouse.get_pressed()[0]:
                pos = pygame.mouse.get_pos()
                if cells[pos[0]//cell_width, pos[1]//cell_width] == 1:
                    cells[pos[0]//cell_width, pos[1]//cell_width] = 0
                else :  
                    cells[pos[0]//cell_width, pos[1]//cell_width] = 1
                update_screen(screen,cells,cell_width)
                pygame.display.update()

        # Dessiner vos cellules ici (remplacez cette partie avec votre propre logique)
        screen.fill(colorGrid)
        if running:
            cells = update(screen,cells,cell_width)
            pygame.display.flip()
        time.sleep(0.001)

        
        clock.tick(fps_max) #The loop won't go faster then 60times/min -> controls the fps


if __name__ == '__main__' :
    args = sys.argv[1:]
    print(args)
    cell_w = int(args[0])
    sys_dim1 = int(args[1])
    sys_dim2 = int(args[2])
    fps = int(args[3])

    main(cell_w, sys_dim1, sys_dim2, fps)