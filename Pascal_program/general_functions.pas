unit general_functions;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, Dialogs,
  general_define_units;

function randomPoisson(mean: real): integer;
procedure ArrayToNegOne(var arr: array of integer);
function ReadMap(mapName: string): Array2Dinteger;
procedure ShowErrorAndExit(const errMsg: string);


implementation

function randomPoisson(mean: real): integer;
  {pseudorandom Poisson distributed number genetrator, Donald Knuth's algorithm}
const
  RESOLUTION = 1000;
var
  k: integer;
  b, l: real;
begin
  //assert(mean > 0, 'mean < 1');
  k := 0;
  b := 1;
  l := exp(-mean);
  while b > l do
  begin
    k := k + 1;
    b := b * random(RESOLUTION) / RESOLUTION;
  end;
  if mean <= 0 then randomPoisson := 0
  else
    randomPoisson := k - 1;
end;

procedure ArrayToNegOne(var arr: array of integer);
var
  i: integer;
begin
  for i := Low(arr) to High(arr) do
    arr[i] := -1;
end;

function ReadMap(mapName: string): Array2Dinteger;
var
  ix, iy, Value: integer;
  mapX, mapY: integer;
begin

  {Read in map}
  Assign(filename, mapName);
  reset(filename);

  readln(filename, mapX, mapY);

  if MapDimX = 0 then
  begin
    MapDimX := mapX;
    MapDimY := mapY;
  end
  else
  if not ((MapDimX = mapX) or (MapDimY = mapY)) then
    ShowErrorAndExit('The current map file has different dimensions from the previous one(s). File: ' + mapName);

  SetLength(Result, Mapdimx + 1, Mapdimy + 1);

  for iy := 1 to Mapdimy do
  begin
    begin
      for ix := 1 to Mapdimx do
      begin
        Read(filename, Value);
        // HabitatMap (and the others) are 'byte' types which is less memory
        // intensive than an integer, but that does mean it can only deal with 0:255
        // There are values in the map of -9999 that represent the sea. As its the same as a barrier, here set to 0!
        if Value < 0 then Result[ix, iy] := 0
        else
          Result[ix, iy] := Value;
      end;
    end;
    readln(filename);
  end;
  Close(filename);

end;

procedure ShowErrorAndExit(const errMsg: string);
begin
  {$IFDEF LCL}
  // If in GUI mode, use MessageDlg for graphical error display
  MessageDlg('Error: ' + errMsg, mtError, [mbOK], 0);
  {$ELSE}
  // In non-GUI mode, use console output
  Writeln('Error: ', errMsg);
  Readln;
  {$ENDIF}

  Halt(1);  // Exit the program with a non-zero code to indicate error
end;


end.

