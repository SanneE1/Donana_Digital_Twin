unit rabbit_input_output_functions;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils,
  general_functions, general_define_units;

procedure ReadClim(breedname, dryname: string);
procedure ReadParameters_Rabbit(paramname: string);
procedure GetAbundanceMap(total_pop: MapOfLists; filename: string; current_sim, year, month: integer);
procedure WritePopsizeMap(filename: string; var Poplist: MapOfLists; dimx, dimy: integer);
procedure WriteArrayToFile(const FileName: string; const arr: array of Integer);
function GetTotalSize(const population_list: MapOfLists): Integer;

implementation

procedure ReadClim(breedname, dryname: string);
var
 C: TextFile;
 Line, Parts: string;
 SplitParts: array of string;
 i, value, RowY, ColX: integer;

begin

SetLength(BreedingMonths_xy, Mapdimx + 1, Mapdimy + 1, (max_years*12));
SetLength(DryMonths_xy, Mapdimx + 1, Mapdimy + 1, (max_years*12));

//Assign(filename, breedname);
//reset(filename);

 WriteLn('Reading Breeding months file');

 AssignFile(C, breedname);
 Reset(C);
 
 while not EOF(C) do
  begin
    ReadLn(C, Line);

    SplitParts := Line.Split([' ']);
    if SplitParts[0] = 'col' then Continue;
    
    ColX := StrToInt(SplitParts[0]);
    RowY := StrToInt(SplitParts[1]);
    
    if (ColX > Mapdimx) or (RowY > Mapdimy) then
      ShowErrorAndExit('There are more columns/rows in the breeding month file than in the map. Check your parameter file to make sure the map files and the climate files match in spatial scale');
    
    SetLength(BreedingMonths_xy[ColX,RowY], (max_years*12) + 1);
    
    for i := 2 to ((max_years*12) + 1) do
    begin
      if i > High(SplitParts) then Break;
      value := StrToInt(SplitParts[i]);
      BreedingMonths_xy[ColX][RowY][i - 1] := value;
    end;
  end;
  
  CloseFile(C);
  
 WriteLn('Reading consecutive dry months file');

 AssignFile(C, dryname);
 Reset(C);
 
 while not EOF(C) do
  begin
    ReadLn(C, Line);

    SplitParts := Line.Split([' ']);
    if SplitParts[0] = 'col' then Continue;
    
    ColX := StrToInt(SplitParts[0]);
    RowY := StrToInt(SplitParts[1]);
    
    if (ColX > Mapdimx) or (RowY > Mapdimy) then
      ShowErrorAndExit('There are more columns/rows in the Consecutive dry months file than in the map. Check your parameter file to make sure the map files and the climate files match in spatial scale');
    
    SetLength(DryMonths_xy[ColX,RowY], (max_years*12) + 1);
    
    for i := 2 to ((max_years*12) + 1) do
    begin
      if i > High(SplitParts) then Break;
      value := StrToInt(SplitParts[i]);
      DryMonths_xy[ColX][RowY][i - 1] := value;
    end;
  end;
  
  CloseFile(C);

end;

procedure ReadParameters_Rabbit(paramname: string);
var
  par_seq: array[1..24] of string;
  val_seq: array of real;
  r, spacePos: integer;
  a, param: string;
  value: real;
begin
  {This function is probably much longer than it needs to be. I just need to make absolutely sure
  that if I at some point change or mess with the param file, I get a warning here, so
  I don't accedentily work with parameter values in the wrong variable!}

   par_seq[1]:= 'max_age';
   par_seq[2]:= 'juv_age';
   par_seq[3]:= 'sub_adult_age';
   par_seq[4]:= 'adult_age';
   par_seq[5]:= 'disp_age';
   par_seq[6]:= 'MeanLitSize';
   par_seq[7]:= 'SdLitSize';
   par_seq[8]:= 'r_int';
   par_seq[9]:= 'r_dens_effect';
   par_seq[10]:= 'r_second_effect';
   par_seq[11]:= 'r_later_effect';
   par_seq[12]:= 'kCapacity_high';
   par_seq[13]:= 'kCapacity_low';
   par_seq[14]:= 'MortP_at_month_old';
   par_seq[15]:= 's_int';
   par_seq[16]:= 's_extra_juv';
   par_seq[17]:= 's_dens_effect';
   par_seq[18]:= 's_food_effect';
   par_seq[19]:= 'lambda';
   par_seq[20]:= 'dens_opt';
   par_seq[21]:= 'sigma';
   par_seq[22]:= 'mapname_rabbit';
   par_seq[23]:= 'breedname';
   par_seq[24]:= 'dryname';

   SetLength(val_seq, High(par_seq)+1);

   Assign(filename, paramname);
   reset(filename);

     for r:=1 to High(par_seq) do
     begin
       readln(filename, a);
       // Find the first space to split the string
      spacePos := Pos(' ', a);

      if spacePos > 0 then
      begin
        // Extract parameter name and convert the rest to a real
        param := Copy(a, 1, spacePos - 1);                      // Get parameter name

        if (param = 'mapname_rabbit') then
          mapname_rabbit := Trim(Copy(a, spacePos + 1, Length(a)))
          else
        if (param = 'breedname') then
          breedname := Trim(Copy(a, spacePos + 1, Length(a)))
          else
        if (param = 'dryname') then
          dryname := Trim(Copy(a, spacePos + 1, Length(a)))
          else
        Val(Trim(Copy(a, spacePos + 1, Length(a))), value);     // Convert value part to real - any integers are converted below to correct type

    if (param = par_seq[r]) then
     val_seq[r] := value
     else
     // stop program and get error message that parameter name not expected
     ShowErrorAndExit('One of the parameter names is not as expected. Check parameter file');
     end
      else ShowErrorAndExit('No space found. Check parameter file');
     end;

     R_max_age             := Round(val_seq[1]);
     R_juv_age             := Round(val_seq[2]);
     R_sub_adult_age       := Round(val_seq[3]);
     R_adult_age           := Round(val_seq[4]);
     R_disp_age            := Round(val_seq[5]);
     R_MeanLitSize         := val_seq[6];
     R_SdLitSize           := val_seq[7];
     R_r_int               := val_seq[8];
     R_r_dens_effect       := val_seq[9];
     R_r_second_effect     := val_seq[10];
     R_r_later_effect      := val_seq[11];

     R_MortP_at_month_old  := val_seq[14];
     R_s_int               := val_seq[15];
     R_s_extra_juv         := val_seq[16];
     R_s_dens_effect       := val_seq[17];
     R_s_food_effect       := val_seq[18];

     if not mcmc_run then
     begin
       kCapacity_high := val_seq[12];
       kCapacity_low  := val_seq[13];

       R_lambda   := val_seq[19];
       R_dens_opt := val_seq[20];
       R_sigma    := val_seq[21];
     end;

     mapname_rabbit := ExpandFileName(mapname_rabbit);
     breedname := ExpandFileName(breedname);
     dryname := ExpandFileName(dryname);

end;

procedure GetAbundanceMap(total_pop: MapOfLists; filename: string; current_sim, year, month: integer);
var
  csvFile: TextFile;
  temp_pop: TList;
  a, x, y: integer;
begin
  AssignFile(csvFile, filename);

  if (current_sim = 1) and (year = 1) and (month = 1) then
    begin
    Rewrite(csvFile);
    // Write header
    WriteLn(csvFile, 'Simulation,Year,Month,N_individuals');
    end
  else append(csvFile);

  a := 0;

  // Write data for each individual
  for x := 1 to MapdimX do
  for y := 1 to MapdimY do
  begin
    temp_pop := RabbitPopulationSpatial[x][y];

    if (temp_pop <> nil) then
     if (temp_pop.Count > 0) then
      begin
      if temp_pop.Count > 1 then RabbitMap[x][y] := RabbitMap[x][y] + 1;          // When (re)calculating breeding habbitat for Lynx, cells are selected for BH and this map will be set back to 0 after
      a := a + temp_pop.Count;
      end;

  end;

  WriteLn(csvFile, current_sim, ',', year, ',', month, ',', a);

  CloseFile(csvFile);
end;


procedure WritePopsizeMap(filename: string; var Poplist: MapOfLists; dimx, dimy: integer);
var
  ix, iy, p: integer;
  outfile: Text;
begin
  Assign(outfile, filename);
  rewrite(outfile);

  // Loop over the arrayData and write each element to the CSV
  for iy := 1 to dimy do
  begin
    for ix := 1 to dimx do
    begin
     if Poplist[ix, iy] <> nil then p := Poplist[ix,iy].Count
     else p := 0;
      // Write each value, followed by a comma, except for the last value in the row
      if ix < dimx then
        Write(outfile, p, ',')
      else
        Write(outfile, p);  // No comma at the end of the row
    end;
    writeln(outfile);  // Move to the next line in the CSV file
  end;

  Close(outfile);
end;


function GetTotalSize(const population_list: MapOfLists): Integer;
var
  i, j: Integer;
  TotalCount: Integer;
begin
  TotalCount := 0;
  for i := 1 to High(population_list) do
    for j := 1 to High(population_list[i]) do
      if Assigned(population_list[i][j]) then
        TotalCount := TotalCount + population_list[i][j].Count;

  Result := TotalCount;
end;

procedure WriteArrayToFile(const FileName: string; const arr: array of Integer);
var
  FileStream: TFileStream;
  i: Integer;
  Line: string;
begin
  // Create or overwrite the file
  FileStream := TFileStream.Create(FileName, fmCreate);
  try
    for i := 0 to High(arr) do
    begin
      Line := IntToStr(arr[i]) + sLineBreak;  // Convert integer to string with newline
      FileStream.Write(Line[1], Length(Line)); // Write to file efficiently
    end;
  finally
    FileStream.Free; // Ensure the file is closed properly
  end;
end;


end.

