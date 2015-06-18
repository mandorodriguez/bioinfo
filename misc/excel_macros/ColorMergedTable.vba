Option Explicit

Sub ColorCels()
   'Colors cells based on their match to the refbase
   
   'Dim snpColors As Dictionary
   'Set snpColors = New Dictionary
   'snpColors("A") = 1
   'snpColors("C") = 2
   'snpColors("T") = 3
   'snpColors("G") = 4
   
   Dim iRowL As Integer
   Dim refbaseCol As Integer
   Dim qIndexStart As Integer
   Dim qIndexEnd As Integer
   Dim stopper As Integer
   Dim numRows As Integer
   Dim numCols As Integer
   Dim refBase As String
   Dim qBase As String
   
   
   Dim r As Integer
   Dim c As Integer
   
   
   'Set up the count as the number of filled rows in the first column of Sheet1.
   numRows = ActiveSheet.UsedRange.Rows.Count
   numCols = ActiveSheet.UsedRange.Rows(1).Cells.Count
   refbaseCol = 4
   qIndexStart = WorksheetFunction.Match("refbase", ActiveSheet.Rows(1), 0) + 1
   qIndexEnd = WorksheetFunction.Match("gene_name", ActiveSheet.Rows(1), 0) - 1
   
   'This loops through all rows
   For r = 2 To numRows
     
     refBase = ActiveSheet.Cells(r, refbaseCol).Value
     
     'This goes through all of the cells
     For c = qIndexStart To qIndexEnd
     
       qBase = ActiveSheet.Cells(r, c).Value
       
       If refBase = qBase Then
       
         ActiveSheet.Cells(r, c).Interior.ColorIndex = 46
       
       End If
       
     
     Next
   
   Next
   
   
   
   stopper = 100

End Sub
