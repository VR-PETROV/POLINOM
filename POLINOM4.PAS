PROGRAM polin;
  uses crt;
  { $R+}{$I+}{$S+}
   VAR nomera,soch:array[1..40] of byte;
      kol_tosh:byte;
      mas_x,mas_y,lagr,pr_lag:array[0..40] of real;
      tabl:array[0..40,0..40] of real;
      i,j,step_polin,nom_slag,step,dl_wyborky,kolvo_chagov,nom_chag,kod:word;
      werh,nig,s:word;
      a,b,sum,chag,x1,znamen,x,polinom,x_v_step,ch,znam_przvd,przvd,pr_pol:real;
      max_przvd,w,max_w:real;

   FUNCTION ZNAK(X:WORD):INTEGER;
    BEGIN
      if x/2=int(x/2) then znak:=1
		      else znak:=-1;
    END;
  PROCEDURE SOCHET(NACH_POZ:BYTE;DL_SOCHET:BYTE;DL_NOMERA:BYTE;UR:BYTE);
    VAR i,j:byte;
	proizv:real;
     BEGIN
       for i:=nach_poz to dl_nomera-dl_sochet+1 do
	begin
	  soch[ur]:=nomera[i];
	  if dl_sochet>1 then
			    sochet(i+1,dl_sochet-1,dl_nomera,ur+1)
			 else
			    begin
			      proizv:=1;
			      for j:=1 to ur do proizv:=proizv*mas_x[soch[j]];
			      sum:=sum+proizv
			    end;

	end;
     END;

  BEGIN
    { ФОРМИРОВАНИЕ ТАБЛИЦЫ ЗАДАЮЩУЮ ФУНКЦИЮ }
    writeln(' ВВЕДИТЕ КОЛИЧЕСТВО ТОЧЕК ПО КОТОРЫМ ИДЕТ ИНТЕРПОЛЯЦИЯ');
    kol_tosh:=0;
    readln(kol_tosh);
    step_polin:=kol_tosh-1;
    writeln('      ');
    writeln('      ');
    for i:=0 to kol_tosh-1 do
     begin
      write(' X',i,' = ');
      read( mas_x[i]);
      write(' Y',i,' = ');
      read(mas_y[i]);
      writeln('      ');
      writeln('      ');
     end;
    a:=mas_x[0];
    b:=mas_x[kol_tosh-1];
    { НАХОЖДЕНИЕ КОЭФФИЦЕНТОВ ПОЛИНОМА ЛАГРАНЖА
      в цикле перебираются все слагаемые полинома лагранжа }

    writeln('ПОИСК КОЭФФИЦЕНТОВ ПОЛИНОМА ЛАГРАНЖА ');
    for nom_slag:=0 to step_polin do
      begin

	{ построение знаменателя i того слагаемого полинома лагранжа
	 здесь програмируется формула (Xi-X0)*(Xi-X1)*...*(Xi-Xi-1)*(Xi-Xi+1)..}
	znamen:=1;
	for j:=0 to step_polin do if nom_slag<>j then znamen:=znamen*(mas_x[nom_slag]-mas_x[j]);

	{ построение массива разрешенных номеров переменной x в слагаемом
	  nom_slag , в  i том слагаемом не должно встречаться Xi }
	i:=1;
	j:=0;
	while j<=step_polin do
	   if nom_slag=j then
			   j:=j+1
			 else
			   begin
			     nomera[i]:=j;
			     i:=i+1;
			     j:=j+1;
			   end;

        { для степени step_polin в таблицу записать коэффицент
          при этой степени }
        tabl[step_polin,nom_slag]:=mas_y[nom_slag]/znamen;

        { вычислить коэффиценты при других степенях в этом слагаемом
	  и записать их в таблицу }
	for step:=step_polin-1 downto 0  do
	  begin
	    dl_wyborky:=step_polin-step;
	    sum:=0;
	    sochet(1,dl_wyborky,step_polin,1);
	    tabl[step,nom_slag]:=znak(dl_wyborky)*sum*mas_y[nom_slag]/znamen;
	  end;
      end;
    { нахождение коэффицентов при степенях многочлена лагранжа
      в результате суммирования по столбцам массива tab  }
      for step:=0 to step_polin do
	begin
	  sum:=0;
	  for j:=0 to step_polin do sum:=sum+tabl[step,j];
	  lagr[step]:=sum;
	end;
    repeat
     { нахождение производной полученного полинома }
     for i:= step_polin-1 downto 0 do pr_lag[i]:=(i+1)*lagr[i+1];

     { вывод полученного полинома }
      writeln(' ПОЛИНОМ   ');
      if step_polin>0 then write(lagr[step_polin]:4:2,'*X^',step_polin);
      for i:=step_polin-1 downto 1 do
	 if lagr[i]>=0 then write('+',lagr[i]:4:2,'*X^',i)
		       else write(lagr[i]:4:2,'*X^',i);
      if lagr[0]>0 then write('+',lagr[0]:4:2);
      if lagr[0]<0 then  write(lagr[0]:4:2);
      readln;



     { вывод производной полинома }
      writeln('     ');
      writeln('     ');
      writeln(' 1ая ПРОИЗВОДНАЯ ПОЛИНОМА  ');
      if step_polin>0 then write(pr_lag[step_polin-1]:4:2,'*X^',step_polin-1);
      for i:=step_polin-2 downto 1 do
	 if pr_lag[i]>=0 then write('+',pr_lag[i]:4:2,'*X^',i)
		       else write(pr_lag[i]:4:2,'*X^',i);
      if pr_lag[0]>0 then write('+',pr_lag[0]:4:2);
      if pr_lag[0]<0 then  write(pr_lag[0]:4:2);
      readln;

     { ввод шага}
     writeln('');
     writeln('');
     repeat
       writeln(' введите шаг');
       readln(chag);
     until (chag<b-a);
     kolvo_chagov:=round( (b-a)/chag);
     writeln('      TАБЛИЦА ЗНАЧЕНИЙ ПОЛИНОМА И ЕГО ПРОИЗВОДНОЙ     ');

     nom_chag:=0;
     while nom_chag<=kolvo_chagov do
       begin
         if (nom_chag/9)=ROUND(nom_chag/9) THEN
           begin
             readln;
             writeln('  ├------------┼-------------┼--------------┤');
	     writeln('  │     X      │  полином(X) │  полином(Х)  │  ');
             writeln('  ├------------┼-------------┼--------------┤ ');
           end;
         x:=a+nom_chag*chag;
      	 { поиск значения полинома в точке X }
         polinom:=0;
	 for step:=0 to step_polin do
	  begin
	    x_v_step:=1;
	    for i:=1 to step do x_v_step:=x_v_step*x;
	    polinom:=polinom+x_v_step*lagr[step];
	  end;

      	 { поиск значения пpоизводной полинома в точке X }
         pr_pol:=0;
	 for step:=0 to step_polin-1 do
	  begin
	    x_v_step:=1;
	    for i:=1 to step do x_v_step:=x_v_step*x;
	    pr_pol:=pr_pol+x_v_step*pr_lag[step];
	  end;

	 write('  │    ',x:4:3);
	 if x>=0 then write(' ');
	 write(' │     ',polinom:4:3);
	 if polinom>=0 then write(' ');
	 write(' │     ',pr_pol:4:3);
	 if pr_pol>=0  then write(' ');
	 writeln('  │   ');

	 writeln('  ├------------┼--------------┤--------------┤ ');
	 nom_chag:=nom_chag+1;
       end;
       while not(keypressed) do a:=a;
       kod:=ord(readkey);
    until kod=27;
   END.

