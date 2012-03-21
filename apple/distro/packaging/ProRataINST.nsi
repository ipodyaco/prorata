!define PRODUCT_NAME "ProRata"
!define PRODUCT_VERSION "1.0 beta"
!define PRODUCT_PUBLISHER "Oak Ridge National Laboratory"
!define PRODUCT_WEB_SITE "http://www.MSProRata.org"
!define PRODUCT_DIR_REGKEY "Software\Microsoft\Windows\CurrentVersion\App Paths\ProRata.exe"
!define PRODUCT_UNINST_KEY "Software\Microsoft\Windows\CurrentVersion\Uninstall\${PRODUCT_NAME}"
!define PRODUCT_UNINST_ROOT_KEY "HKLM"

!include "MUI.nsh"
!include "AddToPath.nsh"

!define MUI_ABORTWARNING
;!define MUI_ICON "${NSISDIR}\Contrib\Graphics\Icons\modern-install.ico"
;!define MUI_UNICON "${NSISDIR}\Contrib\Graphics\Icons\modern-uninstall.ico"
!define MUI_ICON "ProRata32Install.ico"
!define MUI_UNICON "ProRata32Uninstall.ico"

!insertmacro MUI_PAGE_WELCOME
!define MUI_LICENSEPAGE_RADIOBUTTONS
!insertmacro MUI_PAGE_LICENSE "..\ProRata\license.txt"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
;!define MUI_FINISHPAGE_RUN "$INSTDIR\binary\ProRata.exe"
;!define MUI_FINISHPAGE_SHOWREADME "$INSTDIR\Readme.txt"
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_INSTFILES

!insertmacro MUI_LANGUAGE "English"


Name "${PRODUCT_NAME} ${PRODUCT_VERSION}"
OutFile "ProRata_1.0_Setup.exe"
InstallDir "$PROGRAMFILES\ProRata 1.0"
InstallDirRegKey HKLM "${PRODUCT_DIR_REGKEY}" ""
ShowInstDetails show
ShowUnInstDetails show

Section "MainSection" SEC01
  SetOutPath "$INSTDIR\binary"
  SetOverwrite try
  File "..\ProRata\binary\mingwm10.dll"
  File "..\ProRata\binary\ProRata.exe"
  File "..\ProRata\binary\SicForma.exe"
  File "..\ProRata\binary\PRatio.exe"
  File "..\ProRata\binary\ReAdW.exe"
  File "..\ProRata\binary\wolf.exe"
  CreateDirectory "$SMPROGRAMS\ProRata 1.0"
  CreateShortCut "$SMPROGRAMS\ProRata 1.0\ProRata.lnk" "$INSTDIR\binary\ProRata.exe"
  CreateShortCut "$DESKTOP\ProRata 1.0.lnk" "$INSTDIR\binary\ProRata.exe"
  File "..\ProRata\binary\QtCore4.dll"
  File "..\ProRata\binary\QtGui4.dll"
  File "..\ProRata\binary\QtXml4.dll"
  File "..\ProRata\binary\QtNetwork4.dll"
  File "..\ProRata\binary\qwt5.dll"
  SetOutPath "$INSTDIR\docs"
  File "..\ProRata\docs\ProRata-Manual.pdf"
  SetOutPath "$INSTDIR"
  File "..\ProRata\license.txt"
  File "..\ProRata\Readme.txt"
SectionEnd

Section "Add to path"
  Push $INSTDIR\binary
  Call AddToPath
 
  Push "PATH"
  Push $INSTDIR\binary
  Call AddToEnvVar
 
SectionEnd

Section -AdditionalIcons
  WriteIniStr "$INSTDIR\${PRODUCT_NAME}.url" "InternetShortcut" "URL" "${PRODUCT_WEB_SITE}"
  CreateShortCut "$SMPROGRAMS\ProRata 1.0\Manual.lnk" "$INSTDIR\docs\ProRata-Manual.pdf"
  CreateShortCut "$SMPROGRAMS\ProRata 1.0\Website.lnk" "$INSTDIR\${PRODUCT_NAME}.url"
  CreateShortCut "$SMPROGRAMS\ProRata 1.0\Uninstall.lnk" "$INSTDIR\uninst.exe"
SectionEnd

Section -Post
  WriteUninstaller "$INSTDIR\uninst.exe"
  WriteRegStr HKLM "${PRODUCT_DIR_REGKEY}" "" "$INSTDIR\binary\ProRata.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayName" "$(^Name)"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "UninstallString" "$INSTDIR\uninst.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayIcon" "$INSTDIR\binary\ProRata.exe"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "DisplayVersion" "${PRODUCT_VERSION}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "URLInfoAbout" "${PRODUCT_WEB_SITE}"
  WriteRegStr ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}" "Publisher" "${PRODUCT_PUBLISHER}"
SectionEnd


Function un.onUninstSuccess
  HideWindow
  MessageBox MB_ICONINFORMATION|MB_OK "$(^Name) was successfully removed from your computer."
FunctionEnd

Function un.onInit
  MessageBox MB_ICONQUESTION|MB_YESNO|MB_DEFBUTTON2 "Are you sure you want to completely remove $(^Name) and all of its components?" IDYES +2
  Abort
FunctionEnd

Section Uninstall
  Delete "$INSTDIR\${PRODUCT_NAME}.url"
  Delete "$INSTDIR\uninst.exe"
  Delete "$INSTDIR\Readme.txt"
  Delete "$INSTDIR\license.txt"
  Delete "$INSTDIR\binary\qwt5.dll"
  Delete "$INSTDIR\binary\QtXml4.dll"
  Delete "$INSTDIR\binary\QtNetwork4.dll"
  Delete "$INSTDIR\binary\QtGui4.dll"
  Delete "$INSTDIR\binary\QtCore4.dll"
  Delete "$INSTDIR\binary\ProRata.exe"
  Delete "$INSTDIR\binary\SicForma.exe"
  Delete "$INSTDIR\binary\PRatio.exe"
  Delete "$INSTDIR\binary\ReAdW.exe"
  Delete "$INSTDIR\binary\wolf.exe"
  Delete "$INSTDIR\binary\mingwm10.dll"
  Delete "$INSTDIR\docs\ProRata-Manual.pdf"


  Delete "$SMPROGRAMS\ProRata 1.0\Uninstall.lnk"
  Delete "$SMPROGRAMS\ProRata 1.0\Manual.lnk"
  Delete "$SMPROGRAMS\ProRata 1.0\Website.lnk"
  Delete "$DESKTOP\ProRata 1.0.lnk"
  Delete "$SMPROGRAMS\ProRata 1.0\ProRata.lnk"

  Push $INSTDIR\binary
  Call un.RemoveFromPath
 
  Push "PATH"
  Push $INSTDIR\binary
  Call un.RemoveFromEnvVar

  RMDir "$SMPROGRAMS\ProRata 1.0"
  RMDir "$INSTDIR\binary"
  RMDir "$INSTDIR\docs"
  RMDir "$INSTDIR"

  DeleteRegKey ${PRODUCT_UNINST_ROOT_KEY} "${PRODUCT_UNINST_KEY}"
  DeleteRegKey HKLM "${PRODUCT_DIR_REGKEY}"
  SetAutoClose true
SectionEnd
