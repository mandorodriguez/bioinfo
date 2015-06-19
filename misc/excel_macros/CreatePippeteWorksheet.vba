Option Explicit

'This method creates a new worksheet with one column with all of the names from the selection given

Sub CreatePipetteSheet()


Dim rng As Range
Dim s As String
Set rng = Application.Selection
Dim cRow As Range

Dim strDataWSName As String
strDataWSName = ActiveSheet.Name


Dim newWorkSheet As Worksheet
With ActiveWorkbook
    Set newWorkSheet = .Worksheets.Add(After:=.Sheets(.Sheets.Count))
End With

Dim ir, ic As Integer
Dim inwr As Integer
inwr = 1


For ir = 1 To rng.Rows.Count
    Set cRow = rng.Rows(ir)

    For ic = 1 To rng.Columns.Count

        If Application.Sheets(strDataWSName).Cells(cRow.Row, cRow.Columns(ic).Column).Value <> "" Then
        
           newWorkSheet.Cells(inwr, 1).Value = Application.Sheets(strDataWSName).Cells(cRow.Row, cRow.Columns(ic).Column).Value

            inwr = inwr + 1
         End If
         
    Next


Next
 




End Sub