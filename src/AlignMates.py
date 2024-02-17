from cmu_graphics import *
from PIL import Image
import time


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence

class NeedlemanWunsch:
    """
    Performs sequence alignment using the Needleman-Wunsch algorithm
    """
    def __init__(self, seq1, seq2, match, mismatch, gap):
        self.optimPath = []
        self.seq1 = seq1
        self.seq2 = seq2
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.directions = {}
        self.matrix = []
        self.createMatrix()
        self.calculateAlignment()
        # Index for the current position in the optimal path
        self.pathIndex = 0  
        # Time for each cell to be highlighted
        self.pathDisplaySpeed = 0.5  
        # Time in seconds for each cell to be highlighted
        self.lastPathUpdateTime = 0  
        

    def createMatrix(self):
        rows, cols = len(self.seq1.sequence) + 2, len(self.seq2.sequence) + 2
        self.matrix = [[0] * cols for _ in range(rows)]
        for i in range(2, cols):
            self.matrix[0][i] = self.seq2.sequence[i - 2]
            self.matrix[1][i] = (i - 1) * self.gap
        for j in range(2, rows):
            self.matrix[j][0] = self.seq1.sequence[j - 2]
            self.matrix[j][1] = (j - 1) * self.gap

    def calculateAlignment(self):
        rows = len(self.seq1.sequence) + 2
        cols = len(self.seq2.sequence) + 2
        for row in range(2, rows):
            for col in range(2, cols):
                left = self.matrix[row][col - 1] + self.gap
                up = self.matrix[row - 1][col] + self.gap
                diagonal = self.matrix[row - 1][col - 1] + (self.match if self.matrix[0][col] == self.matrix[row][0] else self.mismatch)
                self.matrix[row][col] = max(left, up, diagonal)
                self.directions[(row, col)] = self.calculateDir(left, up, diagonal)

    def calculateDir(self, left, up, diagonal):
        temp = sorted([(left, 'l'), (up, 'u'), (diagonal, 'd')])
        maxValue = temp[-1][0]
        maxValues = [item[1] for item in temp if item[0] == maxValue]
        return tuple(maxValues) if len(maxValues) > 1 else maxValues[0]

    def getNeighbours(self, row, col):
        neighbors = []
        if (row, col) in self.directions:
            for direction in self.directions[(row, col)]:
                if direction == 'd' and row > 1 and col > 1:
                    neighbors.append(('d', row-1, col-1))
                elif direction == 'u' and row > 1:
                    neighbors.append(('u', row-1, col))
                elif direction == 'l' and col > 1:
                    neighbors.append(('l', row, col-1))
        return neighbors

    def reconstructPath(self, path):
        alignedS1 = ''
        alignedS2 = ''
        for direction, i, j in path:
            if direction == 'd':
                alignedS1 = self.seq1.sequence[-i] + alignedS1 
                alignedS2 = self.seq2.sequence[-j] + alignedS2
            elif direction == 'u':
                alignedS1 = self.seq1.sequence[-i] + alignedS1
                alignedS2 = '-' + alignedS2
            else:
                alignedS1 = '-' + alignedS1
                alignedS2 = self.seq2.sequence[-j] + alignedS2
        return (alignedS1, alignedS2)

    def dfs(self, row, col, path, allPaths):
        if row == 1 and col == 1:
            allPaths.append(path)
            return
        for direction, i, j in self.getNeighbours(row, col):
            self.dfs(i, j, path + [(direction, len(self.seq1.sequence)-i+1, len(self.seq2.sequence)-j+1)], allPaths)

    def tracebackDFS(self):
        allPaths = [] 
        self.dfs(len(self.seq1.sequence) + 1, len(self.seq2.sequence) + 1, [], allPaths)
        allAlignments = []
        for path in allPaths:
            alignment = self.reconstructPath(path)
            allAlignments.append(alignment)
        return allAlignments

    def calcAlignScore(self, alignedS1, alignedS2):
        score = 0
        """Calculates the score for a single alignment.
            zip function for parallel iteration: https://www.geeksforgeeks.org/zip-in-python/
            """
        for c1, c2 in zip(alignedS1, alignedS2):
            if c1 == '-' or c2 == '-':
                score += self.gap  # Gap
            elif c1 == c2:
                score += self.match  # Match
            else:
                score += self.mismatch  # Mismatch
        return score

    def calcOptimAlign(self):
        self.optimPath = []  
        allAlignments = self.tracebackDFS()
        maxScore = None
        optimAligns = []
        for alignedS1, alignedS2 in allAlignments:
            score = self.calcAlignScore(alignedS1, alignedS2)
            if maxScore is None or score > maxScore:
                maxScore = score
                optimAligns = [(alignedS1, alignedS2)]
                self.optimPath = [self.convertAlignToPath(alignedS1, alignedS2)]
            elif score == maxScore:
                optimAligns.append((alignedS1, alignedS2))
                self.optimPath.append(self.convertAlignToPath(alignedS1, alignedS2))
        return optimAligns, maxScore

    
    def convertAlignToPath(self, alignedS1, alignedS2):
        path = []
        i, j = len(self.seq1.sequence), len(self.seq2.sequence)   
        for a1, a2 in zip(reversed(alignedS1), reversed(alignedS2)):
            if a1 != '-' and a2 != '-':
                path.append((i + 1, j + 1))  
                i -= 1
                j -= 1
            if a1 == '-' and a2!='-':
                path.append((i+1, j+1))
                j -= 1
            if a2 == '-' and a1!='-':
                path.append((i+1, j+1))
                i -= 1
        return path

    
    def showOptimal(self, app):
        
        optimalAlignments, score = self.calcOptimAlign()
        if score != 0:
            drawLabel(f'Optimal Alignment(s) with score: {score}', 
                    app.width- app.width/3,
                    app.height - app.height/3 - 30, 
                    align='left')
            for i, alignment in enumerate(optimalAlignments):
                alignment0 = alignment[0].replace('-', ' - ')
                alignment1 = alignment[1].replace('-', ' - ')
                drawLabel(f'{alignment0}', 
                        app.width- app.width/3, 
                        app.height - app.height/3 + (6*i+1)*10,
                        align='left')
                drawLabel(f'{alignment1}', 
                        app.width- app.width/3, 
                        app.height - app.height/3 + (6*i+1)*10 +15,
                        align='left')
                


class InfoBubble:
    def __init__(self, app, left, top, width, height, row, col):
        self.app = app
        self.left = left
        self.top = top
        self.width = width
        self.height = height
        self.row = row
        self.col = col
        self.calcDetails = self.calcDetails()

    def draw(self):
        
        cellWidth = 150 
        cellHeight = 40

        
        for i in range(2):
            for j in range(2):
                index = i * 2 + j
                label, value, detail = self.calcDetails[index]
                self.drawCell(self.left + j * cellWidth, self.top + i * cellHeight, cellWidth, cellHeight, label, value, detail)


    def drawCell(self, x, y, width, height, label, value, detail):
       
        drawRect(x, y, width, height, fill="white", border="black")
        textSize = 9  

        labelValText = f"{label}: {value}"
        maxChars = 25   
        detailLines = [detail[i:i+maxChars] for i in range(0, len(detail), maxChars)]

        # Calculate vertical spacing
        totalLines = 1 + len(detailLines)  # 1 for label+value, rest for detail
        lineHeight = height / (totalLines + 1)  # +2 for padding

        drawLabel(labelValText, x + width / 2, y + lineHeight, align="center", size=textSize)
        for i, line in enumerate(detailLines):
            yText = y + (2 + i) * lineHeight  # 2 to account for padding and label+value line
            drawLabel(line, x + width / 2, yText, align="center", size=textSize)


    def calcDetails(self):
        # Calculate each score and the details of the calculation
        details = [("", 0, "") for _ in range(4)] 
        if self.row > 1 and self.col > 1:
            seq1Char = self.app.alignment.seq1.sequence[self.row - 2]
            seq2Char = self.app.alignment.seq2.sequence[self.col - 2]
            diagScore = self.app.alignment.matrix[self.row - 1][self.col - 1]
            matchOrMismatch = self.app.alignment.match if seq1Char == seq2Char else self.app.alignment.mismatch
            diagDetail = f"{diagScore} + {matchOrMismatch} ('match' if {seq1Char} == {seq2Char} else 'mismatch')"
            details[0] = ("Score from Diagonal cell ", diagScore + matchOrMismatch, diagDetail)

        if self.row > 1:
            upScore = self.app.alignment.matrix[self.row - 1][self.col]
            upDetail = f"{upScore} + {self.app.alignment.gap} (Gap)"
            details[1] = ("Score from Upper cell ", upScore + self.app.alignment.gap, upDetail)

        if self.col > 1:
            leftScore = self.app.alignment.matrix[self.row][self.col - 1]
            leftDetail = f"{leftScore} + {self.app.alignment.gap} (Gap)"
            details[2] = ("Score from Left cell ", leftScore + self.app.alignment.gap, leftDetail)

        currentScore = self.app.alignment.matrix[self.row][self.col]
        currentDetail = f"Max of Diagonal, Up, Left"
        details[3] = ("Current cell score", currentScore, currentDetail)

        return details
    

class InputFieldManager:
    def __init__(self, app, runAlignmentFunction):
        self.app = app
        self.runAlignmentFunction = runAlignmentFunction
        self.fields = {
            'sequence1': '',
            'sequence2': '',
            'match': '',
            'mismatch': '',
            'gap': ''
        }
        self.activeField = 'sequence1'
        self.fieldPositions = {
            'sequence1': (app.width - app.width/4, 40),
            'sequence2': (app.width - app.width/4, 90),
            'match': (app.width - app.width/4, 140),
            'mismatch': (app.width - app.width/4, 190),
            'gap': (app.width - app.width/4, 240)
        }
    
    def draw(self, app):
       
        drawLabel('Sequence 1: ', app.width - app.width/4, 50, align='right', size=15)
        drawLabel('Sequence 2: ', app.width - app.width/4, 100, align='right', size=15)
        drawLabel('Match Score: ', app.width - app.width/4, 150, align='right', size=15)
        drawLabel('Mismatch Score: ', app.width - app.width/4, 200, align='right', size=15)
        drawLabel('Gap Score: ', app.width - app.width/4, 250, align='right', size=15)
        
       
        for field in self.fieldPositions:
            if field in ['match', 'mismatch', 'gap']:
                x, y = self.fieldPositions[field]
                drawRect(x, y, 25, 20, fill = None, border='black')
                drawLabel(self.fields[field], x + 5, y+8,align='left', size=15)

            else:
                x, y = self.fieldPositions[field]
                drawRect(x, y, 200, 20, fill = None, border='black')
                drawLabel(self.fields[field], x + 5, y+8,align='left', size=15)
        
    def handleKeyPress(self, key):
        if key == 'tab':  # Switch between input fields
            keys = list(self.fields.keys())
            currentIndex = keys.index(self.activeField)
            self.activeField = keys[(currentIndex + 1) % len(keys)]
        elif key == 'backspace':  # Remove the last character
            self.fields[self.activeField] = self.fields[self.activeField][:-1]
        elif key == 'enter':  # Run Needleman-Wunsch algorithm
            values = self.getValues()
            
            self.runAlignmentFunction(self.app, values['sequence1'], values['sequence2'], values['match'], values['mismatch'], values['gap'])
           
            
        else:  
            if len(self.fields[self.activeField]) < 20:
                self.fields[self.activeField] += key

        
    def getValues(self):
       
        return {
            'sequence1': self.fields['sequence1'],
            'sequence2': self.fields['sequence2'],
            'match': int(self.fields['match']) if not self.fields['match'].isalpha() else "Invalid",
            'mismatch': int(self.fields['mismatch']) if not self.fields['mismatch'].isalpha() else "Invalid",
            'gap': int(self.fields['gap']) if not self.fields['gap'].isalpha() else "Invalid"
        }


def onMouseMove(app, mouseX, mouseY):
    # Activate swipe functionality only if the entire welcome text is displayed
    if app.screen == 'escape' and app.textIndex >= len(app.welcomeText):
        if app.isSwiping:
            if mouseX - app.startMouseX > app.swipeThreshold:
                app.screen = 'right'
                app.isSwiping = False
                app.startMouseX = None
                app.startMouseY = None
            elif abs(mouseY - app.startMouseY) > app.swipeThreshold:
                # Cancel swipe if the vertical movement is too large
                app.isSwiping = False
                app.startMouseX = None
                app.startMouseY = None
        else:
            if app.startMouseX is None and app.startMouseY is None:
                app.startMouseX = mouseX
                app.startMouseY = mouseY
                app.isSwiping = True
    else:
        # Reset swipe detection variables if not on escape screen (i.e. Welcome screen)
        app.isSwiping = False
        app.startMouseX = None
        app.startMouseY = None


    row, col = findWhichCell(app, mouseX, mouseY)
    app.hoveredCell = (row, col) 
    if row is not None and col is not None:
        app.infoBubble = InfoBubble(app, mouseX, mouseY, 100, 100, row + 2, col + 2)
    else:
        app.infoBubble = None

        
def runAlignment(app, seq1, seq2, match, mismatch, gap):
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    
    match = abs(match)

    mismatch = -abs(mismatch)
    gap = -abs(gap)

    app.s1 = Sequence(seq1)
    app.s2 = Sequence(seq2)
    app.alignment = NeedlemanWunsch(app.s1, app.s2, match, mismatch, gap)
    app.needsRedraw = True 



def onKeyPress(app, key):
    # app.inputFieldManager.handleKeyPress(key)
    if key in ['escape', 'space','right']:  
        app.screen = key
      # Switch screen based on key
    else:
        # Existing keypress handling
        app.inputFieldManager.handleKeyPress(key)


# Main application logic
def onAppStart(app):
    app.startMouseX = None
    app.startMouseY = None
    app.isSwiping = False
    app.swipeThreshold = 100
    
    app.bg=Image.open("images/dna.png")
    app.bg=CMUImage(app.bg)
    app.screen = 'escape'
    defaultLength = 8
    app.s1 = Sequence('' * defaultLength)
    app.s2 = Sequence('' * defaultLength)
    app.alignment = NeedlemanWunsch(app.s1, app.s2, 1, -1, -2)
    app.pathColors = ['green', 'blue', 'yellow']
    

    app.inputFieldManager = InputFieldManager(app, runAlignment)
    app.infoBubble = None
    app.boardLeft = app.width * 0.01
    app.boardTop = app.height * 0.3
    app.boardWidth = app.width * 0.55
    app.boardHeight = app.height * 0.55
    app.cellBorderWidth = 2
    app.needsRedraw = False
    app.hoveredCell = None
    app.welcomeText = "Welcome to AlignMates!"
    app.displayedText = ""  
    app.textIndex = 0       # Current index of the text being displayed
    app.textDisplaySpeed = 0.09  # Time in seconds for each character to appear
    app.lastTextUpdateTime = 0  # Time since the last character was added
    
def onStep(app):
    currentTime = time.time()
    if app.screen == 'escape' and app.textIndex < len(app.welcomeText):
        if currentTime - app.lastTextUpdateTime >= app.textDisplaySpeed:
            app.displayedText += app.welcomeText[app.textIndex]
            app.textIndex += 1
            app.lastTextUpdateTime = currentTime
            app.needsRedraw = True

     # Incrementally highlight cells for the traceback path i.e. (slowMo)
    if app.alignment.optimPath and app.alignment.pathIndex < len(app.alignment.optimPath[0]):
        if currentTime - app.alignment.lastPathUpdateTime >= app.alignment.pathDisplaySpeed:
            app.alignment.pathIndex += 1
            app.alignment.lastPathUpdateTime = currentTime
            app.needsRedraw = True



def findWhichCell(app, mouseX, mouseY):
    cellWidth, cellHeight = getCellSize(app)
    row = int((mouseY - app.boardTop) // cellHeight) - 2
    col = int((mouseX - app.boardLeft) // cellWidth) - 2
    if 0 <= row < len(app.s1.sequence) and 0 <= col < len(app.s2.sequence):
        return row, col
    return None, None

def getCellSize(app):
    cellWidth = app.boardWidth / (len(app.s2.sequence) + 2)
    cellHeight = app.boardHeight / (len(app.s1.sequence) + 2)
    return cellWidth, cellHeight

def redrawAll(app):
    if app.screen == 'escape':
        redrawA(app)
    elif app.screen == 'space':
        redrawB(app)
    
    elif app.screen == 'right':
    
        drawLabel("AlignMates",270,100,align="center",size=19,bold=True)
        drawLabel("The dating app for genomic sequences",180,120,align="left",size=12)
        if app.needsRedraw:
            drawBoard(app)
        
        drawBoard(app)
        if app.infoBubble:
            app.infoBubble.draw()

        app.inputFieldManager.draw(app)
        app.alignment.showOptimal(app)
#drawBoard follows similar logic that we implemented in step 1 of tetris in 112
def drawBoard(app):
    for row in range(len(app.s1.sequence) + 2):
        for col in range(len(app.s2.sequence) + 2):
            cellLeft, cellTop = getCellLeftTop(app, row, col)
            cellWidth, cellHeight = getCellSize(app)

            
            fillColor=None
            if app.hoveredCell:
                if app.hoveredCell[0] == row-2 and col == 0:
                    fillColor = "red"
                if app.hoveredCell[1] == col-2 and row == 0:
                    fillColor = "red"

            if app.alignment.optimPath:
                pathSoFar = app.alignment.optimPath[0][:app.alignment.pathIndex]
                if (row, col) in pathSoFar:
                    fillColor = "green" 

            drawRect(cellLeft, cellTop, cellWidth, cellHeight, fill=fillColor, border='black', borderWidth=app.cellBorderWidth)
            drawLabel(app.alignment.matrix[row][col], cellLeft + cellWidth / 2, cellTop + cellHeight / 2)

def redrawA(app):
    drawImage(app.bg, 0, 0, width=app.width, height=app.height)
    drawLabel(app.displayedText, app.width/2, 70, size=50, fill="white")

    # Check if the entire welcome text has been displayed
    if app.textIndex >= len(app.welcomeText):
        drawLabel("Swipe right to find your mate ;)", app.width/2, 480, size=30, fill="white", italic=True)


def redrawB(app):
    """help.png was designed by me using canva(https://www.canva.com/)
    """
    B=Image.open("images/help.png")
    B=CMUImage(B)
    drawImage(B, 0, 0, width=app.width, height=app.height)


def getCellLeftTop(app, row, col):
    cellWidth, cellHeight = getCellSize(app)
    cellLeft = app.boardLeft + col * cellWidth
    cellTop = app.boardTop + row * cellHeight
    return cellLeft, cellTop

def main():
    runApp(width = 1000, height = 600)

main()