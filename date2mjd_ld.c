long double date2mjd_ld(char* ut, char* mjd)
{
  int a,b,year,month,iday,imjd,i;
  double fday;
  char yr[4],mon[2],day[2],hr[2],min[2],sec[2];

  for(i=0;i<4;i++)
    yr[i]=ut[i];
  for(i=0;i<2;i++)
    {
      mon[i]=ut[i+5];
      day[i]=ut[i+8];
      hr[i]=ut[i+11];
      min[i]=ut[i+14];
      sec[i]=ut[i+17];
    }
  year=atoi(yr);
  month=atoi(mon);
  iday=atoi(day);
  fday=(double)(atoi(hr))/24.0+(double)(atoi(min))/24.0/60.0+(double)(atoi(sec))/24.0/3600.0;

  if (month<3) 
    {
      year--;
      month+=12;
    }

  a=floor(year/100.);
  b=2.-a+floor(a/4.);

  if (year<1582) b=0;
  if (year==1582 && month<10) b=0;
  if (year==1582 && month==10 && (double)iday+fday<=4) b=0;

  imjd=floor(365.25*(year+4716))+floor(30.6001*(month+1))+iday+b-2401525;

  sprintf(mjd,"%i.%15.0lf",imjd,fday*1.0e15);

  return (long double)imjd+fday;
}
