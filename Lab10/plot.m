% Wczytanie danych z pliku
data = load('data.txt');

% Załóżmy, że dane są w dwóch kolumnach: x w pierwszej, y w drugiej
x_data = data(:, 1);
y_data = data(:, 2);

% Parametry Lissajous
w1 = 5; 
w2 = 7; 
A1 = 1; 
A2 = 1;
t = 0:2*pi/300:80*pi/min([w1 w2]);
N = 200;

figure; % Utwórz nowe okno wykresu

for i = 1:N
    % Generowanie punktów Lissajous
    x = A1 * sin(w1 * t + i/100);
    y = A2 * sin(w2 * t);
    
    % Rysowanie danych z pliku jako tło
    plot(x_data, y_data, 'k.'); % czarne punkty
    
    hold on;
    % Rysowanie animacji Lissajous
    plot(x, y, 'Color', [i/N 0 (N-i)/N]);
    
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    hold off;
    pause(0.01);
end

