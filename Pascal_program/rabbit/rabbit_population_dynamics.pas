unit rabbit_population_dynamics;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils,
  rabbit_input_output_functions, rabbit_vital_rates,
  general_define_units;

procedure Startpopulation_rabbit;
procedure Rabbit_monthly_demography(current_year, m: integer);
//procedure RunPopSim;

implementation


procedure Startpopulation_rabbit;
var
  a, x, y: integer;
begin
  SetLength(RabbitPopulationSpatial, Mapdimx + 1, Mapdimy + 1);
  SetLength(RabbitMap, MapDimX + 1, MapDimY + 1);

  RabbitPopulation := TList.Create;

  for x := 1 to Mapdimx do
  for y := 1 to Mapdimy do
  begin

  if HabitatMapRabbit[x,y] = 0 then Continue;

  kCapacity := HabitatMapRabbit[x,y];

  RabbitPopulationSpatial[x][y] := TList.Create;

  with RabbitPopulationSpatial[x][y] do
  begin
    for a := 1 to Round(kCapacity * 0.8) do
    begin
      new(Rabbit);

      Rabbit^.age := random(12);  //alternative: Individual^.age:=0;
      if random<0.5 then
        Rabbit^.sex:='f'
      else
        Rabbit^.sex:='m';

      Rabbit^.Pregnant:=False;
      Rabbit^.Lactating:= False;

      RabbitPopulationSpatial[x][y].add(Rabbit);

      //Dispose(Rabbit);
    end;
  end;

  end;
end;

procedure Rabbit_monthly_demography(current_year, m: integer);
var
  b, x, y: integer;
  temp_pop: TList;
begin

  for x := 1 to MapdimX do
  for y := 1 to MapdimY do
  begin

    if HabitatMapRabbit[x,y] = 0 then
    begin
      RabbitPopulationSpatial[x][y] := nil;
      Continue;
    end;

    kCapacity := HabitatMapRabbit[x,y];

    RabbitPopulation := RabbitPopulationSpatial[x][y];
    if RabbitPopulation = nil then Continue;

    RabbitPopulationSize := 0;
    RabbitPopulationSize := RabbitPopulation.Count;

    nbadmonths:= DryMonths_xy[x,y,m];

    if (RabbitPopulationSize > 0) and (tot_m >= RHDData[x,y].NextOutbreakMonth) then
    begin
      ApplyRHD(x, y, tot_m);
      RabbitPopulationSize := RabbitPopulation.Count;
    end;


    Kdensity:=(RabbitPopulationSize)/kCapacity;
        // Zulima's code has (populationsize/5)/kCapacity - but that's because reporting was done at a 5ha scale
        // Currently in the code, we're not calculating # of individuals in the cell, but the Density (#Rabbits/ha)
        // We have a resolution of 25ha per cell (500x500 m).

    reproduction(BreedingMonths_xy[x,y,m]);   //call procedure for reproduction
    RabbitPopulationSize := RabbitPopulation.Count;

    survival;       //call procedure for survival - This needs Kdensity and nbadmonths

        {Aging}
    RabbitPopulationSize := RabbitPopulation.Count;

    if RabbitPopulationSize > 0 then
    with RabbitPopulation do
    for b := 0 to RabbitPopulationSize - 1 do
    begin
      Rabbit := Items[b];
      Rabbit^.Age := Rabbit^.Age + 1;
    end;

    end;
    {With the whole spatial population, start distribution}
    Dispersal;

   //GetAbundanceMap(RabbitPopulationSpatial, output_dir + PathDelim + 'Rabbit_population_size.csv', current_year, m);

  end;

end.

