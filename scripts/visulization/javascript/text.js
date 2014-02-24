function GetColumnId(sTitleLine, sColumn_Name)
{
  var iColumnId = -1;
  sColumns_array = sTitleLine.split("\t");
  for (var i=0; i<sColumns_array.length; i++)
    if (sColumns_array[i] == sColumn_Name)
    {
      iColumnId = i;
      break;
    }
  alert (sTitleLine+ " "+iColumnId.toString());
  return iColumnId;     
}

function parseText (fileText)
{
	var lines = fileText.split(/[\r\n]+/g);
	var currentLine;
	var iLineNum = 0;
	for(var i = 0; i < lines.length; i++)
	{
	  currentLine = lines[i].trim();
	  if (currentLine == "")
	    continue;
	  iLineNum += 1;
    if (iLineNum == 1)
		{
		  var iIdentifier_ColumnId = GetColumnId (currentLine, "Identifier") + 1;
      continue;
    }
    if (iLineNum == 2)
    {
      continue;
    }            
	}
}

window.onload = function() 
{
		var fileInput = document.getElementById('fileInput');
		var fileDisplayArea = document.getElementById('fileDisplayArea');
		fileInput.addEventListener('change', function(e) 
		{
			var file = fileInput.files[0];
			var textType = /text.*/;
			if (file.type.match(textType)) 
			{
				var reader = new FileReader();
        //reader.readAsText(file);
				reader.onload = function(e) 
				{
					//fileDisplayArea.innerText = reader.result;
					//var text = reader.result;

          parseText(reader.result)
				}
				reader.readAsText(file);	
			} else 
			{
				fileDisplayArea.innerText = "File not supported!";
			}
		});
}
