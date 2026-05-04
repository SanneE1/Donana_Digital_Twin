unit rabbit_population_dynamics;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils,
  rabbit_input_output_functions, rabbit_vital_rates,
  general_define_units;

procedure Startpopulation_rabbit;
procedure Rabbit_monthly_demography(a, m: integer);
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

  if HabitatMapRabbit[x,y] = 1 then kCapacity := kCapacity_low;
  if HabitatMapRabbit[x,y] = 2 then kCapacity := kCapacity_high;

  BreedingHabitatMap[x,y] := 1; // initiated first lynx breeding habitat map to start the population

  RabbitPopulationSpatial[x][y] := TList.Create;

  with RabbitPopulationSpatial[x][y] do
  begin
    for a := 1 to Round(kCapacity) * 25 do
    begin
      new(Rabbit);

      Rabbit^.age := random(12);  //alternative: Individual^.age:=0;
      if random<0.5 then Rabbit^.sex:='f' else Rabbit^.sex:='m';
      Rabbit^.Pregnant:=False;

      RabbitPopulationSpatial[x][y].add(Rabbit);

    end;
  end;

  end;
end;

procedure Rabbit_monthly_demography(a, m: integer);
var
  b, x, y: integer;

begin

  for x := 1 to MapdimX do
  for y := 1 to MapdimY do
  begin

    if HabitatMapRabbit[x,y] = 0 then Continue;

    if HabitatMapRabbit[x,y] = 1 then kCapacity := kCapacity_low;
    if HabitatMapRabbit[x,y] = 2 then kCapacity := kCapacity_high;

    RabbitPopulation := RabbitPopulationSpatial[x][y];
    RabbitPopulationSize := RabbitPopulation.Count;

    nbadmonths:= DryMonths_xy[x,y, m+(12*(a-1))];

    Kdensity:=(RabbitPopulationSize)/kCapacity;
        // Zulima's code has (populationsize/5)/kCapacity - but that's because reporting was done at a 5ha scale
        // Currently in the code, by only using populationsize, we're not calculating # of individuals in the cell, but the Density (#Rabbits/ha)
        // We have a resolution of 25ha per cell (500x500 m). If we want to calculate the total number of rabbits, we need to divide popsize by 25, as well as the dens variable in the dispersal
        // Currently however, that makes the model too big to run on my laptop, so I'm first going to see how the calibration goes with #Rabbits/ha
        // Using density instead of abundance, does mean the individual cells are more suceptable to stochasticity

    if BreedingMonths_xy[x,y,m+(12*(a-1))] = 1 then reproduction;   //call procedure for reproduction
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

    GetAbundanceMap(RabbitPopulationSpatial, output_dir + PathDelim + 'Rabbit_population_size.csv', n_sim, a, m);

    if a mod 3 = 0 then 
    WritePopsizeMap((output_dir + PathDelim + 'maps' + PathDelim + 'Rabbit_Population_distribution_' + IntToStr(a) + '.csv'), RabbitPopulationSpatial, mapdimX, mapdimy);

  end;

end.

