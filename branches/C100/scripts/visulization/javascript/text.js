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
  //alert (sTitleLine+ " "+iColumnId.toString());
  return iColumnId;     
}

function parseText (fileText)
{
	var lines = fileText.split(/[\r\n]+/g);
	var currentLine;
	var iLineNum = 0;
    var currentMebInfo_array, currentPeakInfo_array;
    var sCurrentCompound;
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
        iMofZ_ColumnId = GetColumnId (currentLine, "m/z");
        iIntensity_ColumnId = GetColumnId (currentLine, "Intensity");
        iSmiles_ColumnId = GetColumnId (currentLine, "SMILES");
        continue;
      }
      if (currentLine.substring(0, 2)==="M\t") 
      {
          if (iLineNum > 3)
            //do sth for drawing compound and plot
            alert(sCurrentSmiles_array.join(" "));
          currentMebInfo_array = currentLine.split("\t"); 
          sCurrentIdentifier = currentMebInfo_array[iIdentifier_ColumnId];
          var dCurrentMofZ_Identified_array = new Array();  
          var dCurrentIntensity_Identified_array = new Array();
          var dCurrentMofZ_unIdentified_array = new Array();
          var dCurrentIntensity_unIdentified_array = new Array();
          var sCurrentSmiles_array = new Array();
      }else
      {
          currentPeakInfo_array = currentLine.split("\t");
          sCurrentCompound = currentPeakInfo_array[iSmiles_ColumnId];
          if (sCurrentCompound == "NA")
          {
              dCurrentMofZ_unIdentified_array.push(currentPeakInfo_array[iMofZ_ColumnId]);
              dCurrentIntensity_unIdentified_array.push(currentPeakInfo_array[iIntensity_ColumnId]);
          }else
          {
              dCurrentMofZ_Identified_array.push(currentPeakInfo_array[iMofZ_ColumnId]);
              dCurrentIntensity_Identified_array.push(currentPeakInfo_array[iIntensity_ColumnId]);
              sCurrentSmiles_array.push(sCurrentCompound);
          }
      }        
	}
    if (iLineNum > 3)
        //do sth for drawing compound and plot
        alert(sCurrentSmiles_array.join(" "));
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
