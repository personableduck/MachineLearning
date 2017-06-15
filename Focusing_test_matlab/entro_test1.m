image=img_G{1,1};

normalized=image-min(image(:));
normalized=normalized/max(normalized(:));

[y,x]=imhist(normalized);

result = entropy(y)