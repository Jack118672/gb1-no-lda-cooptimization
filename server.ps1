param(
  [switch]$NoOpen
)

$ErrorActionPreference = "Stop"
$root = [System.IO.Path]::GetFullPath($PSScriptRoot)
$listener = $null
$selectedPort = $null

foreach ($port in 8765..8795) {
  try {
    $listener = [System.Net.Sockets.TcpListener]::new([System.Net.IPAddress]::Loopback, $port)
    $listener.Start()
    $selectedPort = $port
    break
  } catch {
    $listener = $null
    continue
  }
}

if (-not $selectedPort) {
  Write-Host "Could not start the local app server on ports 8765-8795."
  Write-Host "Close other local server windows and try again."
  exit 1
}

$url = "http://127.0.0.1:$selectedPort/"
Write-Host "GB1 Multi-Mutant App running at $url"
Write-Host "Press Ctrl+C or close this window to stop the app."
if (-not $NoOpen) {
  Start-Process $url
}

function Get-ContentType($path) {
  $extension = [System.IO.Path]::GetExtension($path).ToLowerInvariant()
  switch ($extension) {
    ".html" { "text/html; charset=utf-8" }
    ".csv" { "text/csv; charset=utf-8" }
    ".json" { "application/json; charset=utf-8" }
    ".js" { "text/javascript; charset=utf-8" }
    ".css" { "text/css; charset=utf-8" }
    ".png" { "image/png" }
    ".jpg" { "image/jpeg" }
    ".jpeg" { "image/jpeg" }
    ".svg" { "image/svg+xml" }
    default { "application/octet-stream" }
  }
}

function Send-Response($stream, $status, $contentType, $bytes) {
  $statusText = if ($status -eq 200) { "OK" } elseif ($status -eq 403) { "Forbidden" } elseif ($status -eq 404) { "Not Found" } else { "Server Error" }
  $header = "HTTP/1.1 $status $statusText`r`nContent-Type: $contentType`r`nContent-Length: $($bytes.Length)`r`nConnection: close`r`nAccess-Control-Allow-Origin: *`r`n`r`n"
  $headerBytes = [System.Text.Encoding]::ASCII.GetBytes($header)
  $stream.Write($headerBytes, 0, $headerBytes.Length)
  $stream.Write($bytes, 0, $bytes.Length)
}

try {
  while ($true) {
    $client = $listener.AcceptTcpClient()
    try {
      $stream = $client.GetStream()
      $reader = [System.IO.StreamReader]::new($stream, [System.Text.Encoding]::ASCII, $false, 1024, $true)
      $requestLine = $reader.ReadLine()
      while (($line = $reader.ReadLine()) -ne $null -and $line.Length -gt 0) { }

      if (-not $requestLine) {
        continue
      }

      $parts = $requestLine.Split(" ")
      $rawPath = if ($parts.Length -ge 2) { $parts[1] } else { "/" }
      $pathOnly = $rawPath.Split("?")[0].TrimStart("/")
      $requestPath = [Uri]::UnescapeDataString($pathOnly)
      if ([string]::IsNullOrWhiteSpace($requestPath)) {
        $requestPath = "index.html"
      }

      $fullPath = [System.IO.Path]::GetFullPath([System.IO.Path]::Combine($root, $requestPath))
      if (-not $fullPath.StartsWith($root, [System.StringComparison]::OrdinalIgnoreCase)) {
        Send-Response $stream 403 "text/plain; charset=utf-8" ([System.Text.Encoding]::UTF8.GetBytes("Forbidden"))
      } elseif (Test-Path -LiteralPath $fullPath -PathType Leaf) {
        Send-Response $stream 200 (Get-ContentType $fullPath) ([System.IO.File]::ReadAllBytes($fullPath))
      } else {
        Send-Response $stream 404 "text/plain; charset=utf-8" ([System.Text.Encoding]::UTF8.GetBytes("Not found"))
      }
    } catch {
      try {
        Send-Response $stream 500 "text/plain; charset=utf-8" ([System.Text.Encoding]::UTF8.GetBytes("Server error: " + $_.Exception.Message))
      } catch { }
    } finally {
      if ($reader) { $reader.Dispose() }
      if ($client) { $client.Close() }
    }
  }
} finally {
  if ($listener) {
    $listener.Stop()
  }
}
